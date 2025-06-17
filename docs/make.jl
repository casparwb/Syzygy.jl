using Documenter
using DocumenterVitepress
# makedocs(;
#     format=DocumenterVitepress.MarkdownVitepress(repo = "...", devbranch = "...", devurl = "dev"),
#     )


# makedocs(sitename="My Documentation")

# makedocs(; 
#     sitename = "Syzygy", 
#     authors = "LuxDL et al.",
#     modules = [DocumenterVitepress],
#     warnonly = true,
#     checkdocs=:all,
#     format=DocumenterVitepress.MarkdownVitepress(
#         repo = https://github.com/casparwb/Syzygy.jl, # this must be the full URL!
#         devbranch = "main",
#         devurl = "dev";
#     ),
#     draft = false,
#     source = "src",
#     build = "build",
#     pages = [
#         "Manual" => [
#             "Get Started" => "manual/get_started.md",
#             "Updating to DocumenterVitepress" => "manual/documenter_to_vitepress_docs_example.md",
#             "Code" => "manual/code_example.md",
#             "Markdown" => "manual/markdown-examples.md",
#             "MIME output" => "manual/mime_examples.md",
#             "DocumenterCitations integration" => "manual/citations.md",
#             "CSS Styling" => "manual/style_css.md",
#             "Authors' badge" => "manual/author_badge.md",
#             "GitHub Icon with Stars" => "manual/repo_stars.md",
#         ],
#         "Developers' documentation" => [
#             "The rendering process" => "devs/render_pipeline.md",
#             "Internal API" => "devs/internal_api.md",
#         ],
#         "api.md",
#     ],
#     plugins = [bib, links],
# )

# DocumenterVitepress.deploydocs(;
#     repo = https://github.com/casparwb/Syzygy.jl, # this must be the full URL!
#     devbranch = "main",
#     push_preview = true,
# )