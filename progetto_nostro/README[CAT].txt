PROGETTO CONTROLLI AUTOMATICI - A.A. 2025/2026
Tipologia B - Traccia 2  : Controllo di un riscaldatore elettrico

GRUPPO  :  [     37     ]

COMPONENTI  : 
- [         Achille Pisani          ]
- [         Alessandro Parmeggiani         ]
- [         Youssef Esam Ebrahim Abou Aiesh          ]


Il presente progetto riguarda la modellazione, l’analisi e il controllo
di un riscaldatore elettrico. Il lavoro è stato svolto in ambiente MATLAB
e Simulink, seguendo le specifiche fornite nella traccia di progetto.

-----------------------------------------------------------------------
REQUISITI SOFTWARE
-----------------------------------------------------------------------
- MATLAB (versione R2022b )
- Control System Toolbox
- Simulink

-----------------------------------------------------------------------
STRUTTURA DEL CODICE
-----------------------------------------------------------------------
Il progetto è articolato in uno script MATLAB principale (che copre l'analisi,
il progetto del controllore e le simulazioni linearizzate) e in un modello
Simulink per la simulazione del sistema non lineare e la verifica finale.

File MATLAB:
- [. . . . . . . . . . . . . . .. ].m

File Simulink punto 4:
- [... . .. . . . . . . . . .].slx

File Simulink punto 5 :
- [... . .. . . . . . . . . .].slx


-----------------------------------------------------------------------
CONTENUTO DELLO SCRIPT MATLAB
-----------------------------------------------------------------------

PUNTO 1 – Modellazione e linearizzazione
- Definizione dei parametri fisici del sistema (masse, calori specifici, ecc.)
- Definizione del modello non lineare (equazioni differenziali)
- Calcolo del punto di equilibrio
- Linearizzazione del modello e calcolo delle matrici di stato (A, B, C, D)
- Analisi degli autovalori e della stabilità in anello aperto

PUNTO 2 – Funzione di trasferimento
- Calcolo della funzione di trasferimento della pianta G(s)
- Analisi dei diagrammi di Bode, poli e zeri
- Definizione delle specifiche frequenziali (bande di disturbo e rumore)

PUNTO 3 – Progetto del regolatore
- Definizione delle specifiche statiche (errore a regime) e dinamiche (sovraelongazione, tempo di assestamento)
- Sintesi del regolatore statico per i requisiti a regime
- Sintesi del regolatore dinamico (polo per attenuazione rumore + rete anticipatrice)
- Analisi delle funzioni di sensitività L(s), S(s), T(s) e verifica dei margini

PUNTO 4 – Test sul modello linearizzato 
- Simulazione in ambiente MATLAB del sistema linearizzato
- Verifica delle risposte a gradino, disturbo sinusoidale e rumore
- Generazione dei grafici comparativi (Scenario lungo e Scenario breve)

{ Punto 4 (simulink) }
-Modello Simulink del sistema linearizzato, utilizzato come supporto
  grafico alle simulazioni MATLAB.

-----------------------------------------------------------------------
PUNTO 5 – Simulazione del sistema non lineare (Simulink)
-----------------------------------------------------------------------
Il modello Simulink implementa il sistema non lineare completo in anello
chiuso, utilizzando il regolatore progettato nel Punto 3. Le simulazioni
consentono di:
- Verificare il comportamento del sistema reale in presenza di ingressi composti
  (riferimento + disturbo + rumore)
- Analizzare la validità locale della linearizzazione
-----------------------------------------------------------------------
OUTPUT GENERATI
-----------------------------------------------------------------------
L’esecuzione del progetto produce:
- Output numerici in Command Window (equilibrio, matrici, margini, errore)
-Figure :
- Figura 1: Diagrammi di Bode della pianta G(s) con vincoli.
- Figura 2: Mappa Poli-Zeri di G(s).
- Figura 3: Diagramma di Bode di L(s) con margini di stabilità evidenziati.
- Figura 4: Bode comparativo L(s), S(s), T(s) con vincoli su disturbo e rumore.
- Figura 5: Bode riassuntivo L(s) e G_e(s).
- Figura 6: Punto 4 - Scenario Lungo (Ingressi w, d a bassa frequenza).
- Figura 7: Punto 4 - Scenario Breve (Ingressi w, d, n ad alta frequenza).
- Figura 8: Punto 4 - Uscite totali (Confronto risposte temporali).