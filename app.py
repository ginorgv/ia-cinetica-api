# ==============================================================================
# API del Compilador Genético TIE - VERSIÓN "CEREBRO" v1.0
# Orden: Obedecer estrictamente. Código completo y sin errores.
# Creado por Raúl Vázquez, implementado bajo orden directa por IA.
# Este backend carga un modelo de conocimiento estructurado ("Cerebro").
# ==============================================================================

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from collections import Counter, defaultdict
import json
import os
import random

# --- 1. INICIALIZACIÓN DE LA APLICACIÓN ---
# La carpeta 'static' contiene el frontend y los modelos de conocimiento.
app = Flask(__name__, static_folder='static')
CORS(app)

# --- 2. DATOS GLOBALES Y MODELOS DE CONOCIMIENTO ---
# El diccionario GENOMIC_BRAINS se llenará al iniciar el servidor.
GENOMIC_BRAINS = {}

CODON_TABLE = { 'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAG', 'AAA'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], '*': ['TAA', 'TAG', 'TGA'] }
CONSTANTES_FISICAS = { 'wG': 0.654, 'wA': 0.481, 'wC': 0.334, 'wT': 0.009, 'K': 60 }

# --- 3. FUNCIONES DEL MOTOR TIE ---

def analyze_kinetics(sequence):
    """Calcula el perfil cinético. Estandariza a formato ADN ('T')."""
    sequence = sequence.upper().replace('U', 'T')
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    profile = []
    for i, codon in enumerate(codons):
        counts = Counter(codon)
        pausa = (counts.get('A', 0) * CONSTANTES_FISICAS['wA'] + counts.get('C', 0) * CONSTANTES_FISICAS['wC'] + counts.get('G', 0) * CONSTANTES_FISICAS['wG'] + counts.get('T', 0) * CONSTANTES_FISICAS['wT']) * 100 + CONSTANTES_FISICAS['K']
        profile.append({'codon_index': i, 'codon': codon, 'pause_ms': round(pausa, 2)})
    return profile

def analyze_grammar_with_brain(sequence, organism):
    """Analiza una secuencia usando el modelo de conocimiento estructurado ("Cerebro")."""
    sequence = sequence.upper().replace('U', 'T')
    brain = GENOMIC_BRAINS.get(organism)
    if not brain:
        return {'error': f'El "Cerebro" para el organismo {organism} no está cargado.'}

    # Extraemos el conocimiento del cerebro
    structural_rules = brain.get("structural_rules", {})
    categories = brain.get("categories", {})
    
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    anomalies = []
    inferred_structure = []

    # Analizamos la secuencia codón por codón
    for i, codon in enumerate(codons):
        # Inferencia de categoría
        category = "DESCONOCIDO" # Default
        if codon in categories.get("INICIADOR", []): category = "INICIADOR"
        elif codon in categories.get("TERMINADOR", []): category = "TERMINADOR"
        elif codon in categories.get("CONTINUADOR", []): category = "CONTINUADOR"
        
        inferred_structure.append(f"{codon}({category})")

        # Detección de anomalías en transiciones (violaciones de la sintaxis)
        if i < len(codons) - 1:
            next_codon = codons[i+1]
            allowed_transitions = structural_rules.get(codon, [])
            if next_codon not in allowed_transitions:
                anomalies.append({
                    'position': i,
                    'subroutine': f"{codon} -> {next_codon}",
                    'type': 'Transición Estructural Anómala'
                })

    return {
        'inferred_structure': " ".join(inferred_structure),
        'anomalies': anomalies
    }

def optimize_sequence(protein_sequence, organism):
    """Diseña una secuencia de ADN usando el conocimiento del "Cerebro"."""
    brain = GENOMIC_BRAINS.get(organism)
    if not brain:
        return "Error: Cerebro para este organismo no disponible."
        
    grammar = brain.get("full_grammar", {})
    codon_freq = defaultdict(int)
    for sub, count in grammar.items():
        for codon in sub.split(" "): codon_freq[codon] += count

    optimized_dna = ""
    for aa in protein_sequence.upper():
        possible_codons = CODON_TABLE.get(aa, [])
        if not possible_codons: continue
        # Elige el codón sinónimo con la frecuencia más alta en la gramática aprendida
        best_codon = max(possible_codons, key=lambda c: codon_freq.get(c, 0))
        optimized_dna += best_codon
        
    return optimized_dna

# --- 4. ENDPOINTS DE LA API ---

@app.route('/')
def serve_index():
    """Sirve la aplicación web principal."""
    return send_from_directory(app.static_folder, 'index.html')

@app.route('/api/analyzer', methods=['POST'])
def handle_analyzer():
    """Endpoint para el Analizador de Gramática y Cinética."""
    try:
        data = request.get_json()
        if not data or 'sequence' not in data: return jsonify({'error': 'No se proporcionó secuencia'}), 400
        sequence = data.get('sequence', '')
        organism = data.get('organism', 'ecoli')
        
        # El informe gramatical ahora es el análisis del "Cerebro"
        grammar_report = analyze_grammar_with_brain(sequence, organism)
        
        return jsonify({
            'kinetic_profile': analyze_kinetics(sequence),
            'grammar_report': grammar_report
        })
    except Exception as e:
        print(f"ERROR en /api/analyzer: {e}")
        return jsonify({'error': f'Error interno del servidor: {e}'}), 500

@app.route('/api/optimizer', methods=['POST'])
def handle_optimizer():
    """Endpoint para el Optimizador de Secuencias."""
    try:
        data = request.get_json()
        if not data or 'protein_sequence' not in data: return jsonify({'error': 'No se proporcionó secuencia de proteína'}), 400
        protein_sequence = data.get('protein_sequence', '')
        organism = data.get('organism', 'ecoli')
        return jsonify({'optimized_dna': optimize_sequence(protein_sequence, organism)})
    except Exception as e:
        print(f"ERROR en /api/optimizer: {e}")
        return jsonify({'error': f'Error interno del servidor: {e}'}), 500

# --- 5. CARGA INICIAL DEL SERVIDOR ---

def load_brains():
    """Carga los modelos de conocimiento ("Cerebros") al iniciar."""
    print("🧠 Cargando 'Cerebros' pre-compilados...")
    try:
        # La ruta se construye de forma segura desde la raíz de la aplicación
        ecoli_path = os.path.join(app.root_path, 'static', 'ecoli_brain.json')
        yeast_path = os.path.join(app.root_path, 'static', 'yeast_brain.json')
        
        with open(ecoli_path, 'r') as f:
            GENOMIC_BRAINS['ecoli'] = json.load(f)
        print(f"✅ 'Cerebro' de E. coli cargado.")
        
        with open(yeast_path, 'r') as f:
            GENOMIC_BRAINS['yeast'] = json.load(f)
        print(f"✅ 'Cerebro' de Levadura cargado.")
            
    except Exception as e:
        print(f"🚨 ERROR CRÍTICO al cargar conocimiento: {e}")
        print("  Asegúrate de que los archivos 'ecoli_brain.json' y 'yeast_brain.json' están en la carpeta 'static'.")

# Ejecutamos la carga de conocimiento al definir la aplicación
load_brains()

if __name__ == '__main__':
    # Esta sección es para pruebas locales. Gunicorn/Render usan su propio método para iniciar.
    app.run(host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))