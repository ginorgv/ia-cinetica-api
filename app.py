# ==========================================================
# Archivo: app.py
# El servidor de tu IA, listo para producción en Render.
# ==========================================================
from flask import Flask, request, jsonify
from flask_cors import CORS
import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
import io

# --- 1. Inicializar la aplicación Flask ---
app = Flask(__name__)
# CORS(app) permite que tu futura página web llame a esta API
CORS(app)

# --- 2. Entrenar y cargar el modelo UNA SOLA VEZ al iniciar el servidor ---
def load_prediction_engine():
    print("🧠 Entrenando y cargando el motor predictivo para E. coli...")
    
    # Los datos verificados de Hussmann et al. (2021)
    datos_reales = """Codon,Pausa_Experimental,Abundancia_ARNt
UUU,81.1,0.0271
UUC,75.0,0.0384
UUA,135.2,0.0125
UUG,111.9,0.0253
UCU,105.1,0.0152
UCC,92.5,0.0212
UCA,120.2,0.0131
UCG,102.8,0.0127
UAU,88.8,0.0234
UAC,76.7,0.0326
CUU,97.3,0.0227
CUC,94.5,0.0211
CUA,166.4,0.0076
CUG,70.9,0.0716
CCU,111.8,0.0107
CCC,188.7,0.0048
CCA,158.4,0.0079
CCG,103.5,0.0191
CAU,93.4,0.0214
CAC,77.5,0.0326
CAA,113.6,0.0206
CAG,83.0,0.0487
CGU,86.8,0.0230
CGC,78.2,0.0298
CGA,244.6,0.0048
CGG,228.4,0.0057
AUU,97.4,0.0323
AUC,77.4,0.0558
AUA,221.7,0.0058
AUG,101.9,0.0427
ACU,102.0,0.0204
ACC,81.9,0.0423
ACA,125.7,0.0123
ACG,93.5,0.0221
AAU,101.4,0.0289
AAC,82.4,0.0494
AAG,79.5,0.0573
AAA,91.8,0.0520
AGU,112.5,0.0145
AGC,84.9,0.0261
AGA,262.8,0.0039
AGG,271.7,0.0037
GUU,92.6,0.0346
GUC,84.4,0.0336
GUA,103.9,0.0226
GUG,97.3,0.0360
GCU,90.2,0.0365
GCC,78.7,0.0560
GCA,86.6,0.0342
GCG,78.6,0.0465
GAU,89.9,0.0368
GAC,78.4,0.0487
GAA,86.5,0.0478
GAG,77.2,0.0633
GGU,119.3,0.0204
GGC,98.6,0.0401
GGA,181.8,0.0084
GGG,160.7,0.0116
UGG,99.9,0.0175
"""
    # He eliminado las columnas de varianza para simplificar
    df = pd.read_csv(io.StringIO(datos_reales)).dropna().reset_index(drop=True)
    df['Codon'] = df['Codon'].str.replace('T', 'U')
    df['S_A'] = df['Codon'].str.count('A')
    df['S_C'] = df['Codon'].str.count('C')
    df['S_G'] = df['Codon'].str.count('G')
    df['S_U'] = df['Codon'].str.count('U')
    df['Escasez_ARNt'] = 1 / df['Abundancia_ARNt']
    
    X = df[['S_A', 'S_C', 'S_G', 'S_U', 'Escasez_ARNt']]
    y = df['Pausa_Experimental']
    model = LinearRegression().fit(X, y)
    tRNA_map = df.set_index('Codon')['Abundancia_ARNt'].to_dict()
    
    print("✅ Motor para E. coli entrenado y listo.")
    return {'model': model, 'tRNA_map': tRNA_map}

PHYSICAL_LAW_MODEL = load_prediction_engine()

# --- 3. Crear el "Endpoint" de la API ---
# Esta es la URL pública que recibirá las secuencias
@app.route('/api/predict/ecoli', methods=['POST'])
def predict_ecoli_profile():
    data = request.get_json()
    if not data or 'sequence' not in data:
        return jsonify({'error': 'No se proporcionó una secuencia en el campo "sequence"'}), 400
    
    sequence = data['sequence']
    model_object = PHYSICAL_LAW_MODEL
    
    try:
        dna_sequence = sequence.upper().replace('T', 'U')
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3) if i+3 <= len(dna_sequence)]
        results = []
        for i, codon in enumerate(codons):
            abundance = model_object['tRNA_map'].get(codon)
            if abundance is None:
                # Omitir codones de parada o desconocidos
                continue
            
            counts = pd.Series(list(codon)).value_counts()
            descriptor = [counts.get('A', 0), counts.get('C', 0), counts.get('G', 0), counts.get('U', 0)]
            tRNA_scarcity = 1 / abundance
            
            input_vector = np.array(descriptor + [tRNA_scarcity]).reshape(1, -1)
            # Pasamos los nombres de las columnas para evitar el UserWarning
            input_df = pd.DataFrame(input_vector, columns=['S_A', 'S_C', 'S_G', 'S_U', 'Escasez_ARNt'])
            score = model_object['model'].predict(input_df)[0]
            
            results.append({'codon_index': i, 'codon': codon, 'predicted_pause_ms': round(score, 2)})
            
        return jsonify(results)
    except Exception as e:
        # Registrar el error real en el servidor para depuración
        print(f"ERROR: {e}")
        return jsonify({'error': 'Ocurrió un error interno al procesar la secuencia.'}), 500