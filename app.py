# ==============================================================================
# API del Compilador Genético TIE (Backend para el MVP) - VERSIÓN FINAL
# Creado por Raúl Vázquez, implementado por IA.
# Este backend carga modelos de gramática pre-compilados para máxima velocidad.
# ==============================================================================

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from collections import Counter, defaultdict
import json
import os
import random

# --- 1. Inicialización de la Aplicación ---
# Servimos los archivos estáticos (index.html, .json) desde la carpeta 'static'
app = Flask(__name__, static_folder='static')
CORS(app)

# --- 2. Datos y Modelos Globales ---
GENOMIC_GRAMMARS = {}
CODON_TABLE = { 'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAG', 'AAA'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], '*': ['TAA', 'TAG', 'TGA'] }
CONSTANTES_FISICAS = { 'wG': 0.654, 'wA': 0.481, 'wC': 0.334, 'wT': 0.009, 'K': 60 }

# --- 3. Funciones del Motor de la IA ---

def analyze_kinetics(sequence):
    sequence = sequence.upper().replace('U', 'T')
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    profile = []
    for i, codon in enumerate(codons):
        counts = Counter(codon)
        pausa = (counts.get('A', 0) * CONSTANTES_FISICAS['wA'] + counts.get('C', 0) * CONSTANTES_FISICAS['wC'] + counts.get('G', 0) * CONSTANTES_FISICAS['wG'] + counts.get('T', 0) * CONSTANTES_FISICAS['wT']) * 100 + CONSTANTES_FISICAS['K']
        profile.append({'codon_index': i, 'codon': codon, 'pause_ms': round(pausa, 2)})
    return profile

def analyze_grammar(sequence, organism):
    sequence = sequence.upper().replace('U', 'T')
    grammar = GENOMIC_GRAMMARS.get(organism, {})
    if not grammar: return []
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    anomalies = []
    for length in [2, 3]:
        for i in range(len(codons) - length + 1):
            subroutine = " ".join(codons[i:i+length])
            if subroutine not in grammar:
                anomalies.append({'position': i, 'subroutine': subroutine, 'type': f'Anomalía (L={length})'})
    return anomalies

def optimize_sequence(protein_sequence, organism):
    grammar = GENOMIC_GRAMMARS.get(organism, {})
    if not grammar: return "Error: Gramática no disponible."
    codon_freq = defaultdict(int)
    for sub, count in grammar.items():
        for codon in sub.split(" "): codon_freq[codon] += count
    optimized_dna = ""
    for aa in protein_sequence.upper():
        possible_codons = CODON_TABLE.get(aa)
        if not possible_codons: continue
        best_codon = max(possible_codons, key=lambda c: codon_freq.get(c, 0))
        optimized_dna += best_codon
    return optimized_dna

# --- 4. Endpoints de la API ---

@app.route('/')
def serve_index():
    return send_from_directory(app.static_folder, 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    return send_from_directory(app.static_folder, path)

@app.route('/api/analyzer', methods=['POST'])
def handle_analyzer():
    data = request.get_json()
    if not data or 'sequence' not in data: return jsonify({'error': 'No se proporcionó secuencia'}), 400
    sequence = data.get('sequence', '')
    organism = data.get('organism', 'ecoli')
    return jsonify({
        'kinetic_profile': analyze_kinetics(sequence),
        'grammar_report': analyze_grammar(sequence, organism)
    })

@app.route('/api/optimizer', methods=['POST'])
def handle_optimizer():
    data = request.get_json()
    if not data or 'protein_sequence' not in data: return jsonify({'error': 'No se proporcionó secuencia de proteína'}), 400
    protein_sequence = data.get('protein_sequence', '')
    organism = data.get('organism', 'ecoli')
    return jsonify({'optimized_dna': optimize_sequence(protein_sequence, organism)})

# --- 5. Carga Inicial del Servidor ---

def load_grammars():
    print("🧠 Cargando modelos de gramática locales desde la carpeta 'static'...")
    try:
        # Buscamos los archivos .json dentro de la carpeta 'static' del servidor
        static_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
        
        with open(os.path.join(static_folder, 'ecoli_grammar.json'), 'r') as f:
            GENOMIC_GRAMMARS['ecoli'] = json.load(f)
        print(f"✅ Conocimiento de E. coli cargado.")
        
        with open(os.path.join(static_folder, 'yeast_grammar.json'), 'r') as f:
            GENOMIC_GRAMMARS['yeast'] = json.load(f)
        print(f"✅ Conocimiento de Levadura cargado.")
        
    except Exception as e:
        print(f"🚨 ERROR CRÍTICO al cargar conocimiento: {e}")
        print("Asegúrate de que los archivos 'ecoli_grammar.json' y 'yeast_grammar.json' están en la carpeta 'static'.")

load_grammars()

if __name__ == '__main__':
    # El puerto es manejado por Gunicorn/Render, pero definimos un default para pruebas locales
    app.run(host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))