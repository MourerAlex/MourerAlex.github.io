# Alex Mourer — Personal Site

Minimal science/blog site. Deploys for free on GitHub Pages.

## Structure

```
site/
├── index.html          ← Homepage
├── css/style.css       ← All styles
├── about/index.html    ← About page
└── blog/
    ├── index.html      ← Blog index (filterable)
    └── posts/
        ├── ate-estimation.html   ← Sample post with LaTeX
        └── ...your posts here
```

## Deploy to GitHub Pages (free, ~5 minutes)

1. Create a free account at github.com
2. Create a new repository named `yourusername.github.io`
3. Upload all files in this folder (drag & drop in the GitHub UI)
4. Go to Settings → Pages → Source: main branch → Save
5. Your site is live at `https://yourusername.github.io`

### Custom domain (optional, ~€10/year)
- Buy a domain (Namecheap, Cloudflare Registrar)
- In GitHub Pages settings, add your custom domain
- Add a CNAME record at your registrar pointing to `yourusername.github.io`

## Adding a new blog post

1. Copy `blog/posts/ate-estimation.html`
2. Update: title, date, tag, content
3. Add the post to `blog/index.html` (copy one `<article>` block)
4. Add it to the homepage `index.html` if it's recent

## Math equations

Uses KaTeX (fast, no server needed).
- Inline: `$your equation$`
- Display block: `$$your equation$$`

## Adding your HTML modules

Drop your existing HTML module files into `blog/posts/` and link to them
from `blog/index.html` like any other post.
