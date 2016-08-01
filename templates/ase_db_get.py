from ase.db import connect

db_path = "{{ db_path }}"
select = "{{ select }}"

db = connect(db_path)

images = []
for d in db.select(select):
    atoms = db.get_atoms(d.id)
    del atoms.constraints
    images += [atoms]
