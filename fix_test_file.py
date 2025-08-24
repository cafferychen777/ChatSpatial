import re

# Read the file
with open('test_data_standardization.py', 'r') as f:
    content = f.read()

# Fix string formatting issues
content = re.sub(r'print\("([^"]*)\{([^}]+)\}([^"]*)"?\)', r'print("\1" + str(\2) + "\3")', content)

# Write back
with open('test_data_standardization.py', 'w') as f:
    f.write(content)

print("Fixed formatting issues in test file")