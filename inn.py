def readMecSol(file_name):

    with open(file_name, "r") as arquivo: #read file
        file = arquivo.readlines()

    def is_number(s): #test if a string is a number
        try:
            float(s)
            return True
        except ValueError:
            pass

        try:
            import unicodedata
            unicodedata.numeric(s)
            return True
        except (TypeError, ValueError):
            pass

        return False

    clean = []
    tags = []

    for i in file:
        if "*" in i: #search for block tags
            tags.append(file.index(i))
        clean.append(i.strip())

    cleaner = {}

    for i in range(len(tags)-1): #get block info for that tag
        cleaner[clean[tags[i]].replace("*", "")] = clean[tags[i]+1:tags[i+1]-1]

    for i in cleaner.keys(): #transform numbers from strings
        temp = []
        for j in cleaner[i]:
            j = [float(x) if is_number(x) else x for x in j.split()] #list comprehensions magic
            temp.append(j)
        cleaner[i] = temp
    return cleaner