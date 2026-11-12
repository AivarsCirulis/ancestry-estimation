path = "aadr_classify_populations.py"
with open(path) as f:
    content = f.read()

old1 = '''    if raw_label in KNOWN_OUTLIER_EXTRA or OUTLIER_PATTERN.search(raw_label):
        return "EXCLUDE", None, None, "QC outlier/discovery/relative subset"'''
new1 = '''    # NOTE: match against `bare` (suffix stripped), not raw_label -- a bare
    # trailing "_o" immediately followed by ".DG" (e.g. "YRI_o.DG") does NOT
    # match "_o$" against the full raw_label, since the string actually ends
    # in ".DG". Found via real data: 13 "_o"-suffixed outlier labels slipped
    # through when this checked raw_label instead of bare.
    if raw_label in KNOWN_OUTLIER_EXTRA or OUTLIER_PATTERN.search(bare):
        return "EXCLUDE", None, None, "QC outlier/discovery/relative subset"'''
assert old1 in content, "old1 not found -- file already patched or differs from expected"
content = content.replace(old1, new1)

old2 = '''    "Turkish.DG": "Europe", "Russian.DG": "Europe",
    "English.DG": "Europe",  # only relevant if included in the British pool'''
new2 = '''    "Turkish.DG": "Europe", "Russian.DG": "Europe",
    "English.DG": "Europe",  # only relevant if included in the British pool
    "IBS_CanaryIslands.DG": "Europe",  # small (n=2) but genuine modern Spanish subgroup, not an outlier'''
assert old2 in content, "old2 not found -- file already patched or differs from expected"
content = content.replace(old2, new2)

with open(path, "w") as f:
    f.write(content)
print("patched OK")
