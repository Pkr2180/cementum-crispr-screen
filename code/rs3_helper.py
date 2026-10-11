"""Rule Set 3 (sequence mode) scoring, patched for lightgbm >= 4.

rs3 0.0.18 ships a LGBMRegressor pickled under lightgbm 3.x; lightgbm 4.x
expects a ``_n_classes`` attribute that the old pickle lacks, so
``model.predict`` raises TypeError. We restore the regressor default (1).
"""
from rs3.seq import load_seq_model, featurize_context

_MODEL = None


def get_model():
    global _MODEL
    if _MODEL is None:
        m = load_seq_model()
        if getattr(m, "_n_classes", None) is None:
            m._n_classes = 1
        _MODEL = m
    return _MODEL


def score_contexts(context_sequences, tracr="Hsu2013", n_jobs=1, chunk=20000):
    """Score 30-mer context sequences (4nt 5' + 20nt protospacer + 3nt PAM + 3nt 3').

    Returns a list of Rule Set 3 sequence scores (higher = more active).
    """
    model = get_model()
    out = []
    for i in range(0, len(context_sequences), chunk):
        feats = featurize_context(context_sequences[i:i + chunk],
                                  sequence_tracr=tracr, n_jobs=n_jobs)
        out.extend(model.predict(feats).tolist())
    return out
