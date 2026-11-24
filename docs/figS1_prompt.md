# Figure S1 Artwork Prompt & Verified Sample-Flow Counts

> **Provenance note**: The final Figure S1 artwork in the manuscript (`docs/figS1_sample_flow.svg`) was generated using Gemini 4 Neon from the prompt and verified counts below (see manuscript LLM statement). All counts are produced programmatically by `scripts/18_sample_flow_audit.R` (and visualized programmatically as a reference vector diagram by `scripts/20_figures.R`).

## Prompt Used to Generate Figure S1

Draw a clean, publication-grade left-to-right sample-flow diagram (174 mm wide, colour-blind-safe Okabe–Ito species palette: AK `#D55E00`, AP `#56B4E9`, BL `#CC79A7`, DI `#009E73`, DV `#0072B2`, LN `#E69F00`, species abbreviations in panels, full italic Latin names in bottom legend, Native `(N)` and Invasive `(I)` marked consistently):

1. **Field Collection (`n = 40` ascidian specimens)**:
   - **Monospecific colonies (`n = 19`)**:
     - `AK` (*Aplidium kottae*, Invasive): 3 (Mount, 29 Apr 2021)
     - `AP` (*Aplidium phortax*, Native): 3 (Mount, 29 Apr 2021)
     - `BL` (*Botrylloides leachii*, Invasive): 3 (Pier A, 15 Mar 2021)
     - `DI` (*Didemnum incanum*, Native): 4 (Mount, 25 Mar 2021)
     - `DV` (*Didemnum vexillum*, Invasive): 3 (Pier A, 15 Mar 2021)
     - `LN` (*Lissoclinum notti*, Native): 3 (Pier A, 15 Mar 2021)
   - **Interaction specimens (`n = 21`; `7` pairings)**:
     - **Native–Invasive (`5` pairings, `n = 15`)**: `DV–DI` 3 (Bridge 18 Feb [2], Mount 25 Mar [1]), `DV–LN` 3 (Pier A 15 Mar), `DV–AP` 3 (Mount 25 Mar), `AK–AP` 3 (Mount 25 Mar), `AK–DI` 3 (Mount 25 Mar).
     - **Invasive–Invasive (`2` pairings, `n = 6`)**: `DV–AK` 3 (Mount 25 Mar), `DV–BL` 3 (Pier A 15 Mar).
   - **Environmental & Controls**: `3` Seawater (Bridge, Pier A, Mount; 11 Aug 2021) + `15` Controls (`13` PCR NTC, `1` extraction blank `S-164`, `1` artificial seawater procedural control `S-163`).

2. **Dissection & Library Prep (`182` total sequenced libraries)**:
   - Monospecific: `19 × 2 = 38` tissue samples (`19` Zooid + `19` Tunic).
   - Interaction: `21 × 6 = 126` tissue samples (`42` Partner A [21 Z + 21 T], `42` Partner B [21 Z + 21 T], `42` Interaction Border [21 Z + 21 T]).
   - Total biological samples: `164` tissue + `3` seawater = **`167` biological libraries** (+ `15` controls = `182` libraries deposited to ENA).

3. **Sequencing, Decontam & Depth QC (`>= 1,000` reads)**:
   - `167` biological libraries (`14` reagent blanks used in plate-batched `decontam`; `92` contaminant ASVs removed, `28` host symbionts rescued [`20` retained post-organelle filter]).
   - **Excluded (`< 1,000` reads, `n = 16`)**:
     - Non-IB (`12`: `7` Z, `5` T): `AP` (`3` Z, `4` T: 4 int, 3 mono), `AK` (`2` Z mono), `BL` (`1` Z mono), `DV` (`1` Z int), `LN` (`1` T int).
     - IB (`4`: `2` Z, `2` T): `AK–AP` (1 T), `DV–AP` (1 Z), `DV–BL` (1 Z), `DV–LN` (1 T).
     - Seawater (`0`).
   - **Retained canonical dataset (`n = 151` libraries, `18,131` ASVs, `6,976,346` reads)**: `73` Zooid + `75` Tunic + `3` Seawater.

4. **Canonical Analysis Sets (`104 + 6 + 38 + 3 = 151`)**:
   - **Primary paired set**: **`52` complete non-IB Zooid–Tunic pairs (`n = 104`)**: `AK` 10, `AP` 5, `BL` 5, `DI` 10, `DV` 17, `LN` 5 (`38` Interaction + `14` Monospecific pairs).
   - **Unpaired non-IB retained**: `n = 6` (`AK` 2 T, `AP` 1 Z, `BL` 1 T, `DV` 1 T, `LN` 1 Z).
   - **Interaction-border (IB) set**: `n = 38` (`19` Zooid + `19` Tunic) $\rightarrow$ **`33` complete IB triplets** (`16` Z + `17` T across `19` encounters: `22` NAT–INV + `11` INV–INV).
   - **Seawater reference**: `n = 3` (Bridge, Pier A, Mount).
