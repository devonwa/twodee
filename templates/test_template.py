#!/usr/bin/env python
###############################################
# test_template.py
# Template for testing Jinja2
# Particularly good to test imports here.

{% if code %}
{{ code }}
{% else %}
print("Hello template world.")
{% endif %}
