from snakemake.utils import validate

_method = config.get("method", "rbfe").strip().lower()

if _method == "abfe":
    validate(config, schema="../../config/schemas/config_abfe.schema.yml")
else:
    validate(config, schema="../../config/schemas/config_rbfe.schema.yml")
