#!/usr/bin/python

# Program to build BK tree
# for faster autocorrection
# @author - Sarthak Sharma <sarthaksharma@gatech.edu>
# Date of Last Modification - 09/07/2018

import distance

class Tree():
    def __init__(self, word):
        self.children = {}
        self.word = word

    def add_child(self, word):
        if self.word:
            ed = distance.hamming(word, self.word)
            if ed not in self.children:
                self.children[ed] = Tree(word)
            else:
                self.children[ed].add_child(word)
        else:
            self.word = word

    def search(self, query, tolerance):
        similar_words = self._search(query, tolerance)
        return similar_words

    def _search(self, query, tolerance):
        results = []
        ed = distance.hamming(self.word, query)
        if ed <= tolerance:
            results.append((ed,self.word))
        for i in range(ed-tolerance,ed+tolerance+1):
            child = self.children.get(i)
            if child is not None:
                results.extend(child._search(query,tolerance))
        return sorted(results)

    def add(self,data):
        if self.data:
            if data < self.data:
                if self.left is None:
                    self.left = Tree(data)
                else:
                    self.left.add(data)
            elif data > self.data:
                if self.right is None:
                    self.right = Tree(data)
                else:
                    self.right.add(data)
        else:
            self.data = data

    def printTree(self):
        if self.children:
            for ed in self.children:
                self.children[ed].printTree()
        print self.word

def main():
    word_list = ["hell","helo","sell"]
    tree = Tree(word_list[0])
    for word in word_list[1:]:
        tree.add_child(word)
    res = tree.search("help",1)
    print res

if __name__ == '__main__':
    main()