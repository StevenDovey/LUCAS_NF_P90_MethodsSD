# DBH model (ingrowth / back-prediction)

For **woody** habits (canopy tree, subcanopy tree, shrub, unknown), missing DBH is filled **per plot** by regressing **response-cycle** DBH on **predictor-cycle** DBH, with a **linear mixed model, random intercept by species**, when the plot has enough species diversity. **One species** or a **failed mixed model** → **plot-only linear regression** (same response and predictor).

Predictions are merged onto the stem record with **dbh1 / dbh2 / dbh3** (third / first / second survey diameters in the usual indexing). **Only missing DBH is replaced**; measured values stay. After each stage, **only missing DBH for that cycle** is filled before the next stage.

**Tree fern, palm, and cabbage tree** use **copy rules** from another cycle, with **non-negative finite checks**, instead of the mixed model in the blocks where those rules run.

## Stage order

Order is **not** “all backward then all forward”.

1. **cycle2←cycle3.** Predictor **dbh3**, response **dbh2**; woody **habit3**, alive **status3**; imputation label **m_2fr3**; then special-habit copy **f_3to2** where coded.

2. **cycle2←cycle1** (forward). Predictor **dbh1**, response **dbh2**; woody **habit1**, alive **status2**; **m_2fr1**; special-habit **f_1to2**.

3. **cycle3←cycle2** (forward). Predictor **dbh2**, response **dbh3**; woody **habit2**, alive **status3**; **m_3fr2**; special-habit **f_2to3** or **f_1to3**.

4. **cycle1←cycle2.** Predictor **dbh2**, response **dbh1**; woody **habit2**, alive **status2**; **pred_dbh1**; special-habit **f_2to1**.

5. **Sparse dead carry:** dead **status2**, **decayclass2** below threshold, **dbh1** missing → carry **dbh2** (and aligned fields), source **sp_21**.

6. **cycle1←cycle3** with **pred_dbh1a** merged with **pred_dbh1**: if one prediction **> 0** and the other not, keep the positive one; if **both > 0**, take the **larger**; both finite but neither **> 0** → combined **0** then non-negative handling; single finite pred → use if **> 0**, else **0** path. Sources **m_1fr2**, **m_1fr3**, **m_1_zero**; then **f_3to1** where coded.

7. **cycle3←cycle1** (forward). Predictor **dbh1**, response **dbh3**; woody **habit1**, alive **status3**; **m_3fr1**.

Because **dbh2** can be imputed **before** step 4, it can act as predictor for **cycle1←cycle2** where applicable.

**Non-negative:** imputed and copied values keep only **finite values ≥ 0**; otherwise **missing**. Not blanket “clamp prediction to zero”; the **only max of two predictions** is in the dual **cycle1** path when **both** are **positive**.

**Diagnostics stack order** for the ingrowth summary: **cycle2←cycle3**, **cycle2←cycle1**, **cycle1←cycle2**, **cycle1←cycle3**, **cycle3←cycle2**, **cycle3←cycle1** (**sp_21** is not a separate regression row there).

**Per-plot regression layer:** mixed model with **random intercept by species**, predict **species-level** adjustment first then **population** line if needed; if **≤1** species or mixed model fails → **plot OLS**; if **<2** complete training rows → no fit; predictions pass through the **same non-negative rule** before merge.
