from snakemake.utils import validate

_method = config.get("method", "rbfe").strip().lower()

if _method == "abfe":
    validate(config, schema="../../config/schemas/config_abfe.schema.yml")
else:
    validate(config, schema="../../config/schemas/config_rbfe.schema.yml")


def _gromacs_cutoff_arg():
    """Cutoff argument for GROMACS min/eq stages, matching RBFE production."""
    if _method != "rbfe":
        return ""
    gromacs_settings = config["production-settings"].get("gromacs-settings", {})
    return f"--cutoff {gromacs_settings.get('cutoff', '9A')}"
