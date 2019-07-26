restart
-- How dictionaries work.

help Dictionary

dictionaryPath

-- symbols get added on the first mutable Dictionary on the path.
select(dictionaryPath, mutable)

-- symbol lookup is done linearly on the dictionaryPath.
dictionary faces
dictionary minimalPrimes

restart
A = new Dictionary
dictionaryPath = prepend(A, dictionaryPath)
A#?"me"
me
A#?"me"
A#"me"
listSymbols A
peek A
dictionary me
dictionary symbol me
me = 42
dictionary "me" -- needs to be a symbol!
str = "me"
value str
dictionaryPath = drop(dictionaryPath, 1)
me
str = "me"
value str
dictionary me
dictionaryPath = prepend(A, dictionaryPath)
me

newPackage "Foo"
dictionaryPath
export "globalsym1"
protect localsym1
globalsym1 = () -> new HashTable from {localsym1 => 42}
peek Foo.Dictionary
peek Foo#"private dictionary"
endPackage "Foo"

dictionary Foo
dictionaryPath

H = globalsym1()
dictionary(first keys H)
debug Foo
dictionary(first keys H)

dictionary symbol me
dictionaryPath
dismiss Foo
dictionaryPath
needsPackage "Foo"

help GlobalDictionary
peek PackageDictionary
peek OutputDictionary
dictionary symbol o42

