import sys

def process_file(to_split):
    x = char_in.split()
    if len(x) > 3:
        print(x[3])
    else:
        print("The file does not contain enough elements.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    char_in = sys.argv[1]
    process_file(char_in)
