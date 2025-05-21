
def print_error(message: str):
    print(f"\x1b[31m[!] Error:\x1b[0m {message}")

def print_info(message: str):
    print(f"\x1b[34m[.] Info:\x1b[0m {message}")

def print_ok(message: str):
    print(f"\x1b[32m[+] Ok:\x1b[0m {message}")

def print_run(message: str):
    print(f"\x1b[33m[%] Run:\x1b[0m {message}")

def print_error_and_exit(message: str):
    print_error(message)
    exit(1)