import json
import os


def load_json(path, default):
    """Ładuje dane z pliku JSON albo zwraca default"""

    if not os.path.exists(path):
        return default

    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except:
        return default


def save_json(path, data):
    """Zapisuje dane do pliku JSON"""

    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    except:
        pass
