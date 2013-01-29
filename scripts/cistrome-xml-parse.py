from xml.etree import ElementTree
import sys
# Usage: python3 cistrome-xml-parse.py ../database/cistrome.xml > ../database/newdb.db
def print_findall(xPath, node,
                  extract=lambda x:x.text.strip(),
                  sep=" ",
                  end="\n"):
    line = ""
    for child in node.findall(xPath):
        line = extract(child) + sep + line

    print(line.strip(),end=end)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("please assign an xml file")
    else:
        tree = ElementTree.parse(sys.argv[1])
        for node in tree.findall('.//motif'):
            print(node.attrib["id"])
            print_findall(".//A", node)
            print_findall(".//G", node)
            print_findall(".//C", node)
            print_findall(".//T", node)
            print()
        # Need to print the last again (to avoid bug)
        print(node.attrib["id"])
        print_findall(".//A", node)
        print_findall(".//G", node)
        print_findall(".//C", node)
        print_findall(".//T", node)
        print()
