# Pushing this repository to GitHub

The folder is complete and ready to publish. Run these commands in PowerShell or Git Bash
from inside this directory.

## 1. Create the repository on GitHub

Either create an empty repository named `cementum-crispr-screen` at
<https://github.com/new> (do **not** add a README, .gitignore or licence — they already
exist here), or use the GitHub CLI in step 3.

## 2. Initialise and commit

```
cd "C:\Users\PRADEEP KUMAR\Documents\Cementum_CRISPR_Screen\cementum-crispr-screen"
git init -b main
git add -A
git commit -m "Cementum CRISPR knockout library: design, off-target verification and manuscript"
```

Set your identity first if git asks:

```
git config --global user.name  "Your Name"
git config --global user.email "you@example.com"
```

## 3. Push

With the GitHub CLI (creates the repository for you):

```
gh repo create cementum-crispr-screen --public --source=. --push
```

Or with a repository you created manually in step 1:

```
git remote add origin https://github.com/<your-username>/cementum-crispr-screen.git
git push -u origin main
```

## Notes

- Total size is about 12 MB across 55 files; no file approaches GitHub's 100 MB limit, so
  Git LFS is not needed.
- `.gitignore` already excludes the large intermediate artefacts (genome index, annotation
  dumps, `.npz` hit tables, Word lock files). Do not remove those entries — the genome
  index and raw hit tables are tens of gigabytes.
- Add a licence before making the repository public. For a resource of this kind, CC BY 4.0
  for the data, figures and manuscript and MIT for the code is a common pairing; GitHub's
  "Add file → Create new file → LICENSE" flow will offer templates.
- If the manuscript is under review, consider pushing as a **private** repository first
  (`--private` with the CLI) and switching to public on acceptance.
