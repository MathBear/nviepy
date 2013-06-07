def getint(string):
    """
    From an input string, gives in output the biggest integer starting 
    from the beginning of the string up to the first non number character.
    Else if the string starts from a character it turns back 1.
    """
    n = 1
    i = 0
    while i<=len(string):
        try: n = int(string[0:i+1])
        except ValueError: return n
        except: print "Unexpected error!"
        i = i + 1
    return n

def countlett(letter,string):
    """Find how many 'letter' are in 'string' with the power convention"""
    r=0
    pos=string.find(letter) # The first 'a' starting from self[0]
    while pos!=-1:
        r+=getint(string[pos+2:])
        pos=string.find(letter,pos+1)
    return r

#=====================================================
def get_substring(st,startpos):
#=====================================================
    """
        Extracts everything between an open and close
        brace from a string.

        INPUT:
            st       -- any string (usually a RootedTree)
            startpos -- an integer such that st[startpos] is
                        the open brace ('{') of the desired substring

        OUTPUT:
            A string containing everything from st[startpos] to
            the corresponding close brace ('}') (inclusive).
    """
    return st[startpos:open_to_close(st,startpos)+1]

#=====================================================
def open_to_close(st,startpos):
#=====================================================
    """ 
        Finds end of a substring enclosed by braces starting at 'startpos'.
        Used by get_substring.
    """

    pos=startpos
    openchar=st[pos]
    if openchar=='{':  closechar='}'
    count=1
    while count>0:
        pos+=1
        if st[pos]==openchar:  count+=1
        if st[pos]==closechar: count-=1
    return pos
