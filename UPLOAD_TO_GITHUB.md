# How to upload this package to GitHub

Upload the **contents of this directory**, including the source code, documentation, recorded R environment, and committed outputs. Do not upload unrelated manuscript build files such as `.aux`, `.log`, or `.synctex.gz` files.

## Recommended repository contents

Keep all of the following:

- `README.md`, which should be the landing page;
- `CODEBOOK.md` and this upload guide;
- both R programs in `code/`;
- `environment/R-session-info.txt`; and
- every file in `outputs/`, including the CSV results and `Sim_fig.pdf`.

Before public release, add the final paper citation and DOI to the top of `README.md` when they are available. Also choose a license in consultation with all coauthors; no license has been imposed by this package.

## Option A: upload using the GitHub website

1. Sign in to GitHub and select **New repository**.
2. Choose a descriptive name, for example `insurance-risk-mrv-section6`.
3. Select public or private visibility as appropriate. Do not ask GitHub to create another README, because this package already contains one.
4. Create the repository, then choose **uploading an existing file**.
5. Drag all contents of `section6-reproducibility/` into the upload area, preserving the `code/`, `environment/`, and `outputs/` directories.
6. Use a commit message such as `Add Section 6 reproducibility code and outputs` and commit the files.
7. Open the repository page and confirm that `README.md` renders and that `outputs/Sim_fig.pdf` and all five CSV files are present.

## Option B: upload using Git in a terminal

Open a terminal and change into `section6-reproducibility/`. Then run the following, replacing the example URL with the empty repository URL supplied by GitHub:

```sh
git init
git add .
git commit -m "Add Section 6 reproducibility code and outputs"
git branch -M main
git remote add origin https://github.com/YOUR-ACCOUNT/YOUR-REPOSITORY.git
git push -u origin main
```

If Git asks for authentication, use GitHub’s documented browser, credential-manager, SSH-key, or personal-access-token workflow; an account password is not accepted for command-line Git operations over HTTPS.

## Final public-release check

From a fresh copy of the uploaded repository, first run the smoke-test command in `README.md`. If adequate computing time is available, then run the full command and compare the regenerated CSV files with the committed outputs. Record any intentional parameter change in both the README and the commit message.
