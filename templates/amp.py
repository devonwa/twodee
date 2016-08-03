#!/usr/bin/env python
###############################################
# amp.py
# Template for building an Amp calculator object

amp_args = {{ amp_args }}

{% if desc_type -%}
# Build descriptor object
from amp.descriptor import {{ desc_type }}
desc = {{ desc_type }}(**{{ desc_args }})
amp_args["descriptor"] = desc
{%- endif %}

{% if reg_type -%}
# Build regression object
from amp.regression import {{ reg_type }}
reg = {{ reg_type }}(**{{ reg_args }})
amp_args["regression"] = reg
{%- endif %}

# Build Amp calculator object
from amp import Amp
calc = Amp(**amp_args)
