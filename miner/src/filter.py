# FIND NEW DATASETS ONLY

def get_new_ids(current_ids, seen_ids):
    return [i for i in current_ids if i not in seen_ids]
