train_args = {{ train_args }}

{% if gs_type %}
# Define global search
from amp import {{ gs_type }}
gs = {{ gs_type }}(**{{ gs_args }})
train_args['global_search'] = gs

{% endif %}
calc.train(images=images,
           **train_args)
