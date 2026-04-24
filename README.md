# Alex Mourer — Personal Site

Static site. No build step. Deployable to GitHub Pages, Cloudflare Pages, Netlify, or any static host.

## Structure
```
index.html                                       Home
about/index.html                                 About
blog/index.html                                  Posts
blog/posts/competing-events/index.html           Series — Competing Events
css/style.css                                    Shared stylesheet
CV_AlexMourer_v5.pdf                             CV
*.png                                            Images
```

## Publish on GitHub Pages
1. Create a repo (e.g. `alexmourer.github.io` or any name).
2. Copy everything in this folder to the repo root.
3. Push. In repo Settings → Pages, set the source to the `main` branch, root.
4. Done. Site goes live at `https://<username>.github.io/` (or `/<repo>/` if not a user-site repo).

## Local preview
```
python3 -m http.server 8000
```
Then open http://localhost:8000

## TODO before going live
- Individual post pages (all `#` links on the home and posts pages).
- Real GitHub repo link in the package section of the Competing Events page.
- Verify LinkedIn / email are current.
