import sys

script = sys.argv[1]

code = open(script).readlines()

found = False
for i, line in enumerate(code):
    if line == "if __name__ == '__main__':\n":
        found = True
        break
if found:
    code = code[:i]
    code.extend(["import unittest.runner\n",
                 "runner = unittest.runner.TextTestRunner()\n",
                 "runner.run(suite())\n"])

script_file = open(script, 'w')
script_file.write(''.join(code))
script_file.close()

