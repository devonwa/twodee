#!/usr/bin/env python
###############################################
# footer.py
# Template for printing standard success message for parsing

{% if print_str %}
# Print data from script.
print({{ print_str }})
{% endif %}

print("=======================================")
print("Successfully completed script.")
print("Exit status: 0")
