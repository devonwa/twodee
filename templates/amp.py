from amp import Amp
from amp.descriptor import {{ desc_type }}
from amp.regression import {{ reg_type }}

desc = {{ desc_type }}(**{{ desc_args }})

reg = {{ reg_type }}(**{{ reg_args }})

calc = Amp(descriptor=desc,
           **{{ amp_args }})
