# Sources and credits

- **invsqrt** — yoyobar's 62-qubit file, byte for byte, from the output of
  https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits (Apache-2.0).

- **dqi** — hengck23's approximate decoder (iteration 014), byte for byte, from
  https://www.kaggle.com/datasets/hengck23/hengck23-dqi-59-57 (Apache-2.0; notebook
  https://www.kaggle.com/code/hengck23/dqi-kmb-at-59-56854), as packaged in Kirill Tymchishin's portfolio
  https://www.kaggle.com/code/kirilltym/quantum-circuits-exact-gauge-5-word-qrom. Decoder lineage as named there:
  SomeJoint's streamed decoder line (https://www.kaggle.com/code/somejoint/quantum-circuits-22-77-local-gm-889-qubit-dqi);
  SomeJoint's research team and GPT-6 Astra Pro; Ahmed Berat Özer; Alex Chilton's 883-qubit Schur decoder, based on
  jday96314's QuantumAI decoder (https://github.com/jday96314/QuantumAI); mrmeows, Ahmed Berat Özer and manwithacat
  (mask projection); Dmitrii Gluzdov, Agniv Banerjee (SomeJoint), ViberTeam, Phuc Thai (DePie).
- **dfthcsos** — Atakan Aldemir's learned angle words on SomeJoint's (Agniv Banerjee's team) R = 8, B = 143, C = 60
  factorization. The fifteen angle circuits, byte for byte, as assembled in Alex Chilton's portfolio (output of
  https://www.kaggle.com/code/alexchilton/22-57-lb-verified-three-task-portfolio, Apache-2.0; per-slot assembly credited
  there to manwithacat): slots 2, 3 and 14 from Dmitrii Gluzdov's
  https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source, slots 0, 1, 4, 5 and 8 ours (from
  our earlier versions), the rest from Atakan Aldemir's
  https://www.kaggle.com/code/atakanaldemir/qcc-888q-repair-and-native-chemistry-co-design (as packaged in
  https://www.kaggle.com/code/kirilltym/quantum-circuits-exact-gauge-5-word-qrom). Loader lineage as named there:
  Atakan Aldemir, heliosli, Dmitrii Gluzdov, on SomeJoint's design; factorization lineage: Michael Rustad, Ahmed Berat
  Özer, Akmal Xodarev (busyaprime), Vladislav Lavrentev (lightningv08), mrmeows. The table weights (W, H2_I, offset_W) are ours, starting from Dmitrii Gluzdov's (version 39 output of
  https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source, Apache-2.0), with
  earlier work starting from his versions 37 and 38, from Sargis Khachatryan's (avazaku; version 1 output
  of https://www.kaggle.com/code/avazaku/free-gap-sqp-889-qubit-dqi, Apache-2.0, built there on our earlier
  version) and from Ahmed Berat Özer's (v14 output of
  https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression, Apache-2.0) and using work from
  sigmaborov's (version 5 output of https://www.kaggle.com/code/sigmaborov/22-85569-lb-exact-gauge-search, Apache-2.0)
  and baidalinadilzhan's (https://www.kaggle.com/code/baidalinadilzhan/quantumai-improved-public-oct29, Apache-2.0;
  its code was used in making ours), all of which refine Atakan's.
- **Organizers** — task definitions, verifiers and kickmix (Google Quantum AI / Qualtran, competition data).
