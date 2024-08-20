import json

def load_config(filename='config.json'):
    try:
        with open(filename, 'r') as f:
            config = json.load(f)
        return config
    except FileNotFoundError:
        print(f"Configuration file {filename} not found.")
        return {}
    except json.JSONDecodeError:
        print(f"Error decoding the configuration file {filename}.")
        return {}


config = load_config()
