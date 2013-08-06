from xml.etree import ElementTree
import sys
# Usage: python3 cistrome-xml-anno-parse.py ../database/cistrome.xml
def print_findall(xPath, node,
                  extract=lambda x:x.text.strip(),
                  sep=" ",
                  end="\n"):

    print(sep.join(map(extract,node.findall(xPath))),end=end)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("please assign an xml file")
    else:
        tree = ElementTree.parse(sys.argv[1])
        for node in tree.findall('.//motif'):
            print(node.attrib["id"],end="\t")
            print_findall(".//symbol",node,end="\t",sep=",")
            print_findall(".//synonym",node,sep=",")
