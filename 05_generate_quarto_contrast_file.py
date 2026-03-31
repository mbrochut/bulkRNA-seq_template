import os
from ruamel.yaml import YAML

# --------------------------------------------------
# Parameters
# --------------------------------------------------
contrast_folder = "./results/contrasts/"
output_folder = "./QUARTO/"
template_file = "template.qmd"

project_title = "HSPC Aging"
author_name = "Maelick Brochut"
split_to_remove = 1
yaml = YAML()
yaml.default_flow_style = False

# --------------------------------------------------
# Setup output folder
# --------------------------------------------------
os.makedirs(output_folder, exist_ok=True)

# --------------------------------------------------
# List contrast files
# --------------------------------------------------
contrast_files = [
    f for f in os.listdir(contrast_folder)
    if f.endswith(".csv")
]

# --------------------------------------------------
# Read QMD template
# --------------------------------------------------
with open(template_file, "r") as f:
    template_content = f.read()

contrast_entries = []

# --------------------------------------------------
# Generate QMD files
# --------------------------------------------------
for contrast_file in contrast_files:
    contrast_name = os.path.splitext(contrast_file)[0]

    split = contrast_name.split("_")
    split = split[split_to_remove:]


    title = " ".join(split)

    contrast_entries.append(
        {"name": contrast_name, "title": title}
    )

    content = template_content
    content = content.replace(
        'contrast = ""',
        f'contrast = "{contrast_name}"',
    )
    content = content.replace(
        'title: ""',
        f'title: "{title}"',
    )
    content = content.replace(
        'author: ""',
        f'author: "{author_name}"',
    )

    output_file = os.path.join(
        output_folder,
        f"{contrast_name}.qmd",
    )

    with open(output_file, "w") as out:
        out.write(content)

    print(f"Generated: {output_file}")

# --------------------------------------------------
# Generate _quarto.yml
# --------------------------------------------------
quarto_config = {
    "project": {
        "type": "website",
        "render": ["QC.qmd"]
        + [f"{e['name']}.qmd" for e in contrast_entries],
    },
    "author": [
        {"name": author_name}
    ],
    "website": {
        "title": project_title,
        "navbar": {
            "left": (
                [{"text": "Quality Control", "href": "QC.qmd"}]
                + [
                    {
                        "text": e["title"],
                        "href": f"{e['name']}.qmd",
                    }
                    for e in contrast_entries
                ]
            )
        },
    },
    "format": {
        "html": {
            "toc": True,
            "toc-location": "left",
        }
    },
}

quarto_yml_file = os.path.join(output_folder, "_quarto.yml")
with open(quarto_yml_file, "w") as f:
    yaml.dump(quarto_config, f)

print(f"Generated: {quarto_yml_file}")
print("All Quarto files generated.")
