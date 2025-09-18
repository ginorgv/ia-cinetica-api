# ==============================================================================
# API del Compilador Genético TIE (Backend para el MVP)
# Creado por Raúl Vázquez, implementado por IA.
# ==============================================================================
from flask import Flask, request, jsonify
from flask_cors import CORS
import pandas as pd
from collections import Counter
import random
import io

# --- 1. Inicialización y Carga de Datos ---
app = Flask(__name__)
CORS(app)

# Diccionario para almacenar las "gramáticas" aprendidas por la IA
GENOMIC_GRAMMARS = {}

# Tabla de codones estándar para el Optimizador
CODON_TABLE = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'], 'K': ['AAG', 'AAA'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], '*': ['TAA', 'TAG', 'TGA']
}

# Constantes físicas para el motor cinético
CONSTANTES_FISICAS = {
    'wG': 0.654, 'wA': 0.481, 'wC': 0.334, 'wT': 0.009, 'K': 60
}

def load_and_learn_grammar(organism_name, url):
    """
    Carga datos y aprende la gramática para un organismo.
    NORMALIZA TODO A FORMATO ARNm (CON 'U').
    """
    print(f"🧠 Aprendiendo gramática para {organism_name}...")
    try:
        df = pd.read_csv(url)
        # --- CORRECCIÓN CLAVE ---
        # Se asegura de que todas las secuencias de entrenamiento usen 'U'.
        sequences = df['sequence'].dropna().str.upper().str.replace('T', 'U').tolist()
        
        counts = Counter()
        for seq in sequences:
            codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
            subroutines_2 = (" ".join(codons[i:i+2]) for i in range(len(codons) - 1))
            subroutines_3 = (" ".join(codons[i:i+3]) for i in range(len(codons) - 2))
            counts.update(subroutines_2)
            counts.update(subroutines_3)
            
        GENOMIC_GRAMMARS[organism_name] = dict(counts)
        print(f"✅ Gramática para {organism_name} aprendida (en formato U).")
    except Exception as e:
        print(f"🚨 ERROR al aprender la gramática para {organism_name}: {e}")

# --- 2. Funciones del Motor de la IA ---

def analyze_kinetics(sequence):
    """Calcula el perfil cinético de una secuencia."""
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    profile = []
    for i, codon in enumerate(codons):
        counts = Counter(codon)
        pausa = (
            counts.get('A', 0) * CONSTANTES_FISICAS['wA'] +
            counts.get('C', 0) * CONSTANTES_FISICAS['wC'] +
            counts.get('G', 0) * CONSTANTES_FISICAS['wG'] +
            counts.get('T', 0) * CONSTANTES_FISICAS['wT']
        ) * 100 + CONSTANTES_FISICAS['K']
        profile.append({'codon_index': i, 'codon': codon, 'pause_ms': round(pausa, 2)})
    return profile

code
Python
# ==========================================================
# REEMPLAZA ESTA FUNCIÓN EN TU app.py
# ==========================================================

def analyze_grammar(sequence, organism):
    """
    Analiza la gramática de una secuencia de entrada del usuario.
    NORMALIZA LA ENTRADA A FORMATO ARNm (CON 'U').
    """
    # --- CORRECCIÓN CLAVE ---
    # Se adapta al usuario. No importa si introduce T o U, se convierte a 'U'.
    sequence = sequence.upper().replace('T', 'U')

    grammar = GENOMIC_GRAMMARS.get(organism, {})
    if not grammar:
        return []

    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    anomalies = []

    # Analizar subrutinas de 2 y 3 codones
    for length in [2, 3]:
        for i in range(len(codons) - length + 1):
            subroutine = " ".join(codons[i:i+length])
            if subroutine not in grammar:
                anomalies.append({
                    'position': i,
                    'subroutine': subroutine,
                    'type': f'Anomalía Gramatical (L={length})'
                })

    # Eliminar redundancias (si [A B C] es raro, no reportar también [A B])
    final_anomalies = []
    positions_reported = set()
    for anomaly in sorted(anomalies, key=lambda x: len(x['subroutine']), reverse=True):
        if anomaly['position'] not in positions_reported:
            final_anomalies.append(anomaly)
            # Marcar todas las posiciones cubiertas por esta anomalía
            for i in range(len(anomaly['subroutine'].split(" "))):
                positions_reported.add(anomaly['position'] + i)
                
    return sorted(final_anomalies, key=lambda x: x['position'])```

def optimize_sequence(protein_sequence, organism):
    """Diseña una secuencia de ADN para una proteína."""
    # En este MVP, usaremos una estrategia simple: elegir el codón más frecuente
    # para ese aminoácido en la gramática del organismo.
    grammar = GENOMIC_GRAMMARS.get(organism, {})
    if not grammar: return "Error: Gramática no disponible."

    # Pre-calcular las frecuencias de codones individuales
    codon_freq = defaultdict(int)
    for sub, count in grammar.items():
        for codon in sub.split(" "):
            codon_freq[codon] += count

    optimized_dna = ""
    for aa in protein_sequence.upper():
        possible_codons = CODON_TABLE.get(aa)
        if not possible_codons: continue

        best_codon = possible_codons[0]
        max_freq = -1
        for codon in possible_codons:
            freq = codon_freq.get(codon, 0)
            if freq > max_freq:
                max_freq = freq
                best_codon = codon
        optimized_dna += best_codon
        
    return optimized_dna

# --- 3. Endpoints de la API ---

@app.route('/api/analyzer', methods=['POST'])
def handle_analyzer():
    data = request.get_json()
    sequence = data.get('sequence', '').upper()
    organism = data.get('organism', 'ecoli')
    
    kinetic_profile = analyze_kinetics(sequence)
    grammar_report = analyze_grammar(sequence, organism)
    
    return jsonify({
        'kinetic_profile': kinetic_profile,
        'grammar_report': grammar_report
    })

@app.route('/api/optimizer', methods=['POST'])
def handle_optimizer():
    data = request.get_json()
    protein_sequence = data.get('protein_sequence', '')
    organism = data.get('organism', 'ecoli')
    
    optimized_dna = optimize_sequence(protein_sequence, organism)
    
    return jsonify({'optimized_dna': optimized_dna})

# --- 4. Carga Inicial del Servidor ---
if __name__ == '__main__':
    # Usamos Gists de GitHub con subconjuntos 
    
    BASE_URL = "https://raw.githubusercontent.com/ginorgv/ia-cinetica-api/main/"
    
    ecoli_data_url = f"{BASE_URL}ecoli_subset.csv"
    yeast_data_url = f"{BASE_URL}yeast_subset.csv"
    
    load_and_learn_grammar('ecoli', ecoli_data_url)
    load_and_learn_grammar('yeast', yeast_data_url)
    app.run(host='0.0.0.0', port=10000)
