def levenshtein(ref, seq):
    import numpy as np

    if len(ref) < len(seq):
        return levenshtein(seq, ref)
    if len(seq) == 0:
        return 'ZERO SEQ' #len(ref)

    #print('LEN REF = %i, LEN SEQ = %i' % (len(ref), len(seq)))

    ref = np.array(tuple(ref))
    seq = np.array(tuple(seq))
    #print('REF = %i\nSEQ = %i' % (len(ref), len(seq)))

    previous_row = np.arange(seq.size + 1)
    for s in ref:
        current_row = previous_row + 1
#        print('s = %s\nPREVIOUS ROW = %s\nNEW ROW = %s' %(s, previous_row, current_row))
        current_row[1:] = np.minimum(current_row[1:], np.add(previous_row[:-1], seq != s))
        current_row[1:] = np.minimum(current_row[1:], current_row[0:-1] + 1)
        previous_row = current_row

    return previous_row[-1]

if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2])
   
