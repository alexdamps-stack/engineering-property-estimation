# Entwicklungsplan – Engineering Tools / Property Estimation

> **Historische Fassung.** Dieses Dokument ist die ursprüngliche Planungsgrundlage vom 1. November 2026 und wird nicht mehr fortgeschrieben. Bestandsaufnahme (§2) und Architekturhinweis (§4) beschreiben den Stand vor Einrichtung des Repositorys. Maßgeblich für Umfang, Arbeitspakete und Abnahmekriterien ist der englische Plan [`docs/development-plan.md`](docs/development-plan.md); die Prüfanforderungen stehen in [`docs/testing-strategy.md`](docs/testing-strategy.md).

Stand: 1. November 2026. Neuer Planungsvorschlag auf Grundlage der zugänglichen Projektchats und der lokalen Bestandsaufnahme.

## 1. Ziel und Leitentscheidungen

Eine eigenständige, in Rust entwickelte Anwendung für Stoffdaten, Gemische und stationäre Prozesssimulation aufbauen. Der thermodynamische Kern bleibt unabhängig von Bedienoberfläche und Flowsheet. Zunächst müssen Modelle und Gleichgewichtsberechnungen nachvollziehbar funktionieren; darauf entsteht eine kleine grafische Prozesssimulation.

Aus den jüngsten Nutzerbeiträgen übernommen:

- Die Implementierung soll in Rust erfolgen.
- Excel ist vorerst zurückgestellt.
- Neben dem Rechenkern ist eine Flowsheet-Oberfläche gewünscht.
- Für Kohlenwasserstoffe sollen GERG oder PR-BM berücksichtigt werden.
- Die Trennung zwischen Property Packages und Simulator soll sich am besprochenen DWSIM-Konzept orientieren.
- Spline-basierte Stoffdatenberechnung ist als spätere Beschleunigung zu untersuchen.

Noch zu klären: Im Chat zu Rust wurde auch eine Nutzung in DWSIM erwähnt. Ob eine konkrete DWSIM-Anbindung zusätzlich zur eigenständigen Oberfläche benötigt wird, bleibt ein offener Integrationspunkt. Die interne Rust-API wird zuerst definiert; eine Fremdsprachen-Schnittstelle wird erst nach dieser Entscheidung entworfen.

## 2. Bestandsaufnahme und Aussagegrenzen

| Bereich | Nachweisbarer Stand | Konsequenz |
|---|---|---|
| Lokaler Projektordner | Bei der Prüfung leer; kein Git-Repository | Projektbasis und Versionsverwaltung fehlen lokal |
| Angefügte Artefakte dieses Chats | Keine | Kein direkt verfügbarer Quellcodeanhang |
| Früherer Prototyp | Im verwandten Projektchat als Python-Projekt ThermoCore v0.3.0 beschrieben | Als fachliche Vorarbeit behandeln; keinen bestehenden Rust-Kern voraussetzen |
| Modelle des Prototyps | SRK, PR, PR-BM, NRTL, Flüssigvolumenkorrelationen und GERG-Anbindung über teqp laut früherem Bericht | Formulierungen, Parameter und Schnittstellen vor einer Portierung prüfen |
| Frühere Prüfung | Berichtet wurden 53 Tests und 111 GERG-Vergleiche über 37 Gaszustände | Historischer Bericht, in dieser Sitzung nicht erneut ausgeführt |
| Phasengleichgewicht | Frühere Berichte nennen fehlenden vollständigen Flash und fehlende globale Stabilitätsabsicherung | Im neuen Plan zentrale Entwicklungsarbeit |
| Rust-Solver und Research | Nutzerwünsche sind lesbar; mehrere frühere Assistentenantworten nur als Inhaltsreferenzen verfügbar | Fertigstellungsstand daraus nicht ableiten |

Der frühere Quellcode wurde unter `/workspace/scratch/2a408b01b5b5` erzeugt; im Chat ist `ThermoCore_v0_3_0_Source.zip` verlinkt. Dieser damalige Remote-Pfad ist kein hier verfügbarer lokaler Projektbestand. Die Bestandsaufnahme umfasst die zugänglichen Chattexte und einen ergänzenden Einblick in den zugehörigen früheren Ausführungs-Chat, keine erneute Codeprüfung oder vollständige Prüfung aller Projektdateien in ChatGPT.

## 3. Empfohlener erster Produktumfang

### Erste technische Version

- Stoff- und Gemischdefinitionen mit klarer Mengenbasis, SI-Einheiten und Absolutdruck.
- PR und PR-BM als getrennt auswählbare Modelle.
- Einphasige Zustände: Dichte, Kompressibilitätsfaktor, Fugazitätskoeffizienten und kalorische Größen mit dokumentiertem Referenzzustand.
- TP-Flash für ausgewählte Kohlenwasserstoffgemische einschließlich Stabilitätsprüfung im definierten Testbereich.
- Kommandozeilenprogramm mit reproduzierbaren Ein- und Ausgaben.

### Erste bedienbare Version

- Native Rust-Oberfläche mit Stoffauswahl, Gemischdefinition und Ergebnisansicht.
- Flowsheet mit Feed, Mischer, Heizer/Kühler, Ventil und Zweiphasenabscheider.
- PH-Flash als Voraussetzung für Ventil und energiebasierte Berechnungen.
- Speichern und Laden eines versionierten Projektformats.
- Ein durchgängiges Demonstrationsprojekt: Gasgemisch mischen, kühlen, entspannen und abscheiden.

GERG-2008 ist ein eigener, wichtiger Entwicklungsstrang. Eine erste PR-BM-Demonstration darf entstehen, während GERG implementiert wird. Eine als GERG-fähig bezeichnete Veröffentlichung setzt dagegen den geprüften nativen GERG-Pfad voraus.

Zurückgestellt: Excel, Dynamik, Reaktionen, Elektrolyte, Destillationskolonnen, vollständige Erdölcharakterisierung und allgemeine Mehrflüssigphasensysteme. Kerosin und Diesel bleiben zunächst explizite Surrogatgemische.

## 4. Architekturvorschlag

Die folgenden Namen sind geplante Rust-Crates und noch keine vorhandenen Dateien.

| Modul | Verantwortung |
|---|---|
| `thermo-types` | Einheiten, Stoffidentitäten, Gemische, Zustände, Mengenbasen und Fehler |
| `thermo-data` | Versionierte Stoffparameter, Wechselwirkungsparameter, Quellen und Gültigkeitsbereiche |
| `thermo-models` | PR, PR-BM, später GERG und NRTL; Zustandsgrößen und Ableitungen |
| `thermo-solvers` | Nullstellensuche, Stabilitätsprüfung, TP-/PH-/später PS-Flash und Diagnostik |
| `property-packages` | Zusammenstellung kompatibler Modelle, Daten und Berechnungsfähigkeiten |
| `flowsheet-core` | Ströme, Apparate, Verbindungen, Bilanzierung und Berechnungsreihenfolge |
| `thermo-cli` | Reproduzierbare Beispiele, Stapelberechnung und Regression |
| `thermo-desktop` | Native Oberfläche und Flowsheet-Editor; ruft dieselbe Rechen-API auf |

Ein Property Package legt die Verfahren für Gleichgewicht, kalorische Größen und gegebenenfalls gesonderte Flüssigdichte fest. Zulässige Kombinationen und Referenzzustände müssen geprüft sein. Eine Dichtekorrelation darf nicht stillschweigend thermodynamische Ableitungen oder Fugazitäten eines anderen Modells ersetzen.

Berechnungsergebnisse enthalten Modell- und Datenversion, Phasen, Einheiten, Gültigkeitsstatus und Solverdiagnostik. Fehlende Parameter, nicht unterstützte Phasengebiete und Nichtkonvergenz sind unterscheidbare Ergebnisse. Ein vorheriger erfolgreicher Zustand darf nach einer Eingabeänderung nicht weiter als aktuell erscheinen.

## 5. Arbeitspakete und Abnahmekriterien

| ID | Priorität | Arbeitspaket | Abhängigkeit | Fertig, wenn … |
|---|---|---|---|---|
| A0 | P0 | Frühere Artefakte sichern und inventarisieren | Keine | Quellcode, Formulierungen, Parameterquellen und Tests sind gesichert oder als nicht verfügbar markiert |
| A1 | P0 | Rust-Projektbasis und Rechenvertrag | A0 oder dokumentierter Neuaufbau | Workspace, Formatierung, Linting und CI stehen; SI-, Mengenbasis- und Fehlerkonventionen sind festgelegt |
| A2 | P0 | Stoffdaten und Referenzfälle | A1 | Startkomponenten sind eindeutig identifiziert; Datenherkunft, Referenzzustände und zulässige Bereiche sind erfasst |
| A3 | P0 | PR und PR-BM | A2 | Wurzeln, Fugazitäten und kalorische Größen bestehen unabhängige Vergleiche; BM-Formulierung und Ableitungen sind dokumentiert |
| A4 | P0 | Stabilität und TP-Flash | A3 | Ein-/Zweiphasenfälle, Grenzfälle und Nichtkonvergenz werden korrekt behandelt; definierte Bilanz- und Gleichgewichtsgrenzen werden eingehalten |
| A5 | P0 | PH-Flash und Energiegrundlage | A4 | Enthalpiereferenz und ideale Beiträge sind konsistent; TP→PH-Rückrechnungen und Ventilfälle bestehen |
| A6 | P1 | Nativer GERG-2008-Pfad | A2; Integration nach A4 | Komponenten, Parameter und Ableitungen sind gegen eine unabhängige Referenz geprüft; unterstützte Zustände sind explizit dokumentiert |
| A7 | P1 | NRTL und Flüssigkeitseigenschaften | A2, A4, A5 | Gerichtete Parameter, Temperaturabhängigkeit, Dampfmodell und kalorische Kopplung sind geprüft; fehlende Daten werden sichtbar abgewiesen |
| A8 | P1 | Stationärer Flowsheet-Kern | A5 | Demonstrationsprozess läuft ohne GUI mit geschlossenen Stoff- und Energiebilanzen |
| A9 | P1 | Native Oberfläche und Flowsheet-Editor | API-Prototyp nach A3; Rechenabnahme nach A8 | Nutzer kann das Demonstrationsprojekt erstellen, berechnen, ändern, speichern und laden |
| A10 | P2 | Recycles und größere Netze | A8, A9 | Ein definierter Rückführungsfall konvergiert reproduzierbar; Abbruch und Fehlerpfad sind geprüft |
| A11 | P2 | Tabellen und Splines | Stabile direkte Modelle und Profiling | Beschleunigung und Fehler sind im freigegebenen Bereich gemessen; außerhalb erfolgt direkter Modellaufruf oder klarer Fehler |
| A12 | P2 | Externe Integration | Stabile API und geklärter Zielumfang | Ein konkreter DWSIM- oder später Excel-Anwendungsfall ist mit passender Schnittstelle nachgewiesen |

P0 bildet den kritischen Pfad zur belastbaren technischen Version. P1 ergänzt die bedienbare Anwendung und die vorgesehenen Stoffsysteme. P2 folgt nach einem funktionierenden Gesamtprozess.

## 6. Solverstrategie

1. Skalierte Residuen, physikalische Grenzen, Abbruchbedingungen und reproduzierbare Startwerte definieren.
2. Robuste skalare Nullstellensuche für Dichte- und äußere Temperaturprobleme bereitstellen; Bereiche und Phasenwechsel berücksichtigen.
3. Rachford–Rice für vorgegebene K-Werte gesondert prüfen. Das allein ist noch kein thermodynamischer Flash.
4. Für kubische Modelle TP-Flash mit Fugazitätsgleichheit und gedämpfter sukzessiver Substitution beginnen.
5. TPD-Stabilitätsprüfung mit mehreren Startzusammensetzungen integrieren; nach einer Phasenteilung erneut Stabilität prüfen. Endliche Startversuche nicht als mathematischen Beweis globaler Stabilität ausweisen.
6. Newton-Verfahren mit Schrittbegrenzung beziehungsweise Line Search gezielt ergänzen, sobald Ableitungen und Vergleichsfälle belastbar sind.
7. PH-Flash zunächst als kontrollierte äußere Temperaturiteration über TP-Flash entwickeln; schwierige Phasenübergänge gesondert testen.
8. Kritische Nähe, verschwindende Komponenten und nahezu identische Phasen explizit behandeln. Eine kleine Schrittweite allein ist kein Konvergenznachweis.

Automatische Differentiation, Trust-Region-Verfahren, Homotopie und Reduktionsmethoden werden an konkreten schwierigen Fällen verglichen. Eine universelle Solverarchitektur wird erst bei nachgewiesenem Bedarf ausgebaut.

Die vorgeschlagene Trennung von Stabilitätsanalyse und Phasenteilung ist auch in der aktuellen CoolProp-Dokumentation für Gemische beschrieben: [CoolProp Mixtures](https://coolprop.org/fluid_properties/Mixtures).

## 7. Modell- und Validierungsmatrix

| Stoffsystem | Vorgesehener Pfad | Zentrale Prüfung |
|---|---|---|
| Erdgas | PR-BM für erste Flash-Fälle, GERG-2008 als eigener hochwertiger Modellpfad | Dichte und kalorische Größen; gleiche Zusammensetzung und Referenzzustände; Phasengebiet separat qualifizieren |
| Kohlenwasserstoffgemische | PR-BM | Alpha-Funktion, Wechselwirkungsparameter, VLE und Ableitungen |
| Kerosin-/Dieselsurrogate | PR-BM; gesonderte Flüssigdichte nur mit dokumentierter Kopplung | Modellfehler von Parameter- und Surrogatfehler trennen |
| Reines Methanol | Geprüfter Reinstoffpfad beziehungsweise unabhängige Referenz | Flüssigdichte, Dampfdruck und Enthalpie |
| Methanol/Wasser | NRTL mit ausdrücklich festgelegtem Dampfmodell | Parameterkonventionen, VLE, Temperaturableitungen und kalorische Konsistenz |
| Feuchtes Rauchgas | Idealgas als begrenzter Ausgangsfall; Kondensation gesondert modellieren | Nasse/trockene Zusammensetzung, Wasserabscheidung und Bilanz |

Startdatensatz: Methan, Ethan, Propan, Stickstoff, CO₂, Sauerstoff, Wasser und Methanol; schwere Surrogatkomponenten erst mit geprüften Daten ergänzen.

Für jeden Referenzfall werden Eingabe, Zusammensetzungsbasis, Modellvariante, Parameter, Referenzsoftwareversion, Quellen und erwartete Größen gespeichert. PR wird nur mit passender PR-Referenz verglichen; ein Standard-PR-Ergebnis ist keine unmittelbare PR-BM-Referenz. Modellvergleich, Implementierungsprüfung und experimentelle Validierung werden getrennt ausgewiesen.

Vorgeschlagene anfängliche numerische Abnahmegrenzen für reguläre, gut konditionierte Fälle:

- Zusammensetzungssummen und normierte Komponentenbilanz: absolute Abweichung höchstens `1e-10` beziehungsweise `1e-8`.
- Zweiphasengleichgewicht für aktive Komponenten: maximaler Betrag von `ln(f_i,L/f_i,V)` höchstens `1e-7`.
- Energiebilanz: skalierter Fehler höchstens `1e-6`, mit dokumentierter absoluter Untergrenze für die Skalierung.
- Eigenschaftsvergleiche: pro Größe absolute und relative Toleranzen vor der Implementierung anhand der Referenzqualität festlegen.

Diese Werte sind Entwicklungsziele, keine bereits erreichten Ergebnisse oder universellen Genauigkeitszusagen. Für Spurenkomponenten und kritische Zustände sind eigene Kriterien und Diagnosefälle nötig.

## 8. GERG, Oberfläche und Splines: Entscheidungspunkte

**GERG:** Die frühere teqp-Anbindung erfüllt allein keine vollständig native Rust-Implementierung. teqp und gegebenenfalls weiterer unabhängig geprüfter Referenzcode können außerhalb der ausgelieferten Anwendung Vergleichswerte erzeugen. Vor der nativen Umsetzung sind Formulierungsumfang, idealer und residualer Anteil, Parameterherkunft sowie Nutzungsbedingungen zu erfassen. Die offizielle teqp-Dokumentation verweist auf GERG-2004 und GERG-2008: [GERG-Dokumentation](https://pages.nist.gov/teqp-docs/en/stable/models/GERG.html).

**Oberfläche:** Ein kurzer Rust-GUI-Prototyp prüft Canvas, Verbindungen, Eigenschaftseditor, Hintergrundberechnung und Windows-Paketierung. Die Frameworkentscheidung fällt anhand dieses Prototyps. Das Rechenmodell darf keine Abhängigkeit vom GUI-Framework erhalten. Die Trennung von Simulation und Property Packages entspricht dem dokumentierten DWSIM-Konzept: [DWSIM Simulation](https://dwsim.org/help/main/06-simulation.html).

**Splines:** Zuerst Laufzeitprofile direkter Berechnungen erstellen. Ein erster Tabellenversuch verwendet ein festes Gemisch in einem klar begrenzten einphasigen Bereich. Fehler von Werten und Ableitungen getrennt messen; Phasengrenzen und kritische Bereiche ausnehmen. Interpolierte Zustandsgrößen müssen thermodynamische Konsistenz berücksichtigen. Veränderliche Zusammensetzung, hoher Speicherbedarf und belastbare Rückfallregeln sind eigene Arbeitspunkte. Aus der früheren Dyad-Frage wird keine bereits getroffene technische Entscheidung abgeleitet.

## 9. Zeitliche Orientierung

Unverbindliche Aufwandsschätzung für eine erfahrene Vollzeitentwicklung; nach A0–A3 anhand realer Ergebnisse neu schätzen. Zeiträume sind Entwicklungswochen ab tatsächlichem Start, keine zugesagten Kalendertermine.

| Etappe | Grober Aufwand | Ergebnis |
|---|---|---|
| A0–A2: Bestand, Basis und Daten | 1–2 Wochen | Prüffähiger Ausgangsstand und Referenzfälle |
| A3: PR/PR-BM | 2–4 Wochen | Reproduzierbare einphasige Rechnungen |
| A4–A5: TP-/PH-Flash | 4–8 Wochen | Technische Version im festgelegten Stoff- und Zustandsbereich |
| A8–A9: Flowsheet und Oberfläche | 4–7 Wochen | Bedienbare PR-BM-Demonstration |
| A6: Native GERG-Implementierung | Zusätzlich zunächst 4–10+ Wochen | Separater Aufwand nach Daten- und Formulierungsprüfung |
| A7: NRTL und Flüssigkeitspakete | Zusätzlich zunächst 3–6 Wochen | Erweiterung auf definierte polare Systeme |

Die erste PR-BM-Anwendung liegt damit grob bei 11–21 Entwicklungswochen. GERG, NRTL, kritische Zustände und Integrationen erweitern diesen Umfang. Vorliegende belastbare Artefakte können Aufwand reduzieren; fehlende Parameter und zusätzliche Phasengebiete können ihn erhöhen.

## 10. Unmittelbar nächster Arbeitsblock

1. Den früheren ThermoCore-Stand und verfügbare Solverberichte sichern; fehlende Artefakte im Inventar kennzeichnen.
2. Einen verbindlichen Beispielprozess und dessen Komponenten, Zustandsbereich und benötigte Ergebnisgrößen festhalten.
3. PR-BM-Formulierung, Stoffdaten, Parameter und kalorischen Referenzzustand dokumentieren.
4. Rust-Workspace, Datenformat und CLI-Eingabeformat einrichten.
5. Einen Erdgasfall und einen binären Kohlenwasserstofffall mit unabhängigen Referenzen festlegen und implementieren.
6. Erst nach bestandenem einphasigem Vergleich mit TP-Flash und Stabilitätsprüfung fortfahren.

Erstes Abnahmeergebnis: Eine reproduzierbare Rust-CLI berechnet beide Fälle mit ausgewiesenem Modell, Datenstand, Einheiten und Diagnose und besteht den dokumentierten Referenzvergleich.

## 11. Verwendete Projektquellen

- [Thermodynamic Solver Research](https://chatgpt.com/c/6adc861e-4de0-83ed-b04b-688703c6190d): Rust-Vorgabe und Zurückstellen von Excel; weitere Assistentenantworten teilweise nur als Inhaltsreferenzen verfügbar.
- [Rust statt C# – Projektchat](https://chatgpt.com/c/6ae3c138-1120-83ed-930c-4760a2fc3ccc): eigenständige Anwendung, Property-Package-/Flowsheet-Trennung, gewünschte Oberfläche und Frage zu Spline-Verfahren.
- [Excel Addin Entwicklung](https://chatgpt.com/c/6ad7db8d-3c14-83eb-9b63-e96821e6ebe2): ältere Entwicklungsdiskussion; sichtbare Inhalte erlauben keine belastbare Rekonstruktion der damaligen Umsetzung.
- [Excel Tool Structure Plan](https://chatgpt.com/c/6ad67a9d-dd54-83eb-980f-b2c48b98df6a): Python-Prototyp, Beispielgemische, PR-BM-/GERG-Anforderungen und berichtete Tests.
- Zugehöriger früherer Ausführungs-Chat `01a155d0-d5a2-7574-9f7e-7a2de227c07e`: ergänzende Protokolle zur damaligen Python-Implementierung.
- Lokale Dateisystemprüfung und Prüfung der Chat-Anhänge am 1. November 2026.

Die externen Dokumentationslinks stützen technische Orientierung. Architektur, Prioritäten, Abnahmegrenzen und Aufwand sind Vorschläge dieses neuen Plans.
