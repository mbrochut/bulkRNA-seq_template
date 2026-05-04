import os
import shutil
from ruamel.yaml import YAML
import yaml
# --------------------------------------------------
# Load config
# --------------------------------------------------
with open("config.yaml") as f:
    config = yaml.safe_load(f)

# --------------------------------------------------
# Parameters from config
# --------------------------------------------------
params = config.get("05_generate_quarto", {})
general = config.get("general", {})

contrast_folder = params.get("contrast_folder", "./results/contrasts/")
output_folder = params.get("output_folder", "./QUARTO/")
template_dir = params.get("template_dir", "./Quarto_template/")

project_title = params.get("project_title", "My Project")
author_name = params.get("author_name", "Author")
split_to_remove = params.get("split_to_remove", 1)

organism = general.get("organism", "Human")
# --------------------------------------------------
# LOAD MODULES
# --------------------------------------------------
modules_config = params.get("modules", {})
modules = {}

for key, mod in modules_config.items():
    if not mod.get("include", False):
        continue  # skip disabled modules

    modules[mod["name"]] = mod

yaml = YAML()
yaml.default_flow_style = False

# --------------------------------------------------
# Setup output folder
# --------------------------------------------------
os.makedirs(output_folder, exist_ok=True)

# --------------------------------------------------
# Always copy utils.py
# --------------------------------------------------
if os.path.exists("utils.py"):
    shutil.copy("utils.py", os.path.join(output_folder, "utils.py"))
    print("Copied: utils.py")

# --------------------------------------------------
# Helper functions
# --------------------------------------------------
def build_contrast_title(name, split_to_remove=1):
    parts = name.split("_")[split_to_remove:]
    return " ".join(parts)


def generate_qmd(template_path, output_path, context):
    with open(template_path, "r") as f:
        content = f.read()

    for key, value in context.items():
        content = content.replace(f'{key} = ""', f'{key} = "{value}"')
        content = content.replace(f'{key}: ""', f'{key}: "{value}"')

    with open(output_path, "w") as f:
        f.write(content)


# --------------------------------------------------
# Get contrast files
# --------------------------------------------------
contrast_files = [
    f for f in os.listdir(contrast_folder)
    if f.endswith(".csv")
]

# --------------------------------------------------
# Main generation
# --------------------------------------------------
navbar = []
render_files = []

for module_name, cfg in modules.items():

    if not cfg["include"]:
        continue

    # ------------------------------------------
    # SINGLE FILE MODULE (e.g. QC)
    # ------------------------------------------
    if cfg["type"] == "single":
        src = os.path.join(template_dir, cfg["file"])
        dst = os.path.join(output_folder, cfg["file"])

        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"Copied: {cfg['file']}")

            navbar.append({
                "text": module_name,
                "href": cfg["file"],
            })

            render_files.append(cfg["file"])
        else:
            print(f"⚠️ Missing: {src}")

    # ------------------------------------------
    # MENU MODULE (per contrast)
    # ------------------------------------------
    elif cfg["type"] == "menu":
        template_path = os.path.join(template_dir, cfg["template"])

        if not os.path.exists(template_path):
            print(f"⚠️ Missing template: {template_path}")
            continue

        entries = []

        for contrast_file in contrast_files:
            contrast_name = os.path.splitext(contrast_file)[0]
            title = build_contrast_title(contrast_name, split_to_remove)
            filename = f"{cfg['prefix']}{contrast_name}.qmd"
            output_path = os.path.join(output_folder, filename)

            context = {
                "contrast": contrast_name,
                "title": title,
                "organism": organism,
                "model": cfg['model']
            }

            generate_qmd(template_path, output_path, context)

            entries.append({
                "text": title,
                "href": filename,
            })

            render_files.append(filename)

            print(f"Generated: {filename}")

        if entries:
            navbar.append({
                "text": module_name + f" ({cfg['model']})",
                "menu": entries,
            })

# --------------------------------------------------
# Generate _quarto.yml
# --------------------------------------------------
quarto_config = {
    "project": {
        "type": "website",
        "render": render_files,
    },
    "author": [{"name": author_name}],
    "website": {
        "title": project_title,
        "navbar": {
            "left": navbar
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