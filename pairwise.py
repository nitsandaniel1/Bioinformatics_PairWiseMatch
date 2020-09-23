from copy import copy
import sys
import networkx as nx
from datetime import datetime


class HSP:
    def __init__(self, seq1_start, seq1_end, seq2_start, seq2_end, score):
        self.seq1_start = seq1_start
        self.seq1_end = seq1_end
        self.seq2_start = seq2_start
        self.seq2_end = seq2_end
        self.score = score

    def size(self):
        return self.seq1_end - self.seq1_start

    def diagonal(self):
        return self.seq1_start - self.seq2_start

    def __eq__(self, other):
        """equals function"""
        if isinstance(other, self.__class__):
            is_seq1_eq = self.seq1_start == other.seq1_start and self.seq1_end == other.seq1_end
            is_seq2_eq = self.seq2_start == other.seq2_start and self.seq2_end == other.seq2_end
            return is_seq1_eq and is_seq2_eq and self.score == other.score
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        return hash(str((self.seq1_start, self.seq1_end, self.seq2_start, self.seq2_end, self.score)))

    def __str__(self):
        return (
            f'Score: {self.score}\n'
            f'Seq1: [{self.seq1_start}, {self.seq1_end}]\n'
            f'Seq2: [{self.seq2_start}, {self.seq2_end}]\n'
        )

    def __repr__(self):
        return self.__str__()


k = 10
T = 50
X = 50


def get_alphabet(path):  # checked
    alphabet = ''
    with open(path) as f:
        chars = f.readline().strip().split()
    for char in chars:
        alphabet += char
    return alphabet


def read_scoring_matrix(path):
    scoring_matrix = {}

    with open(path) as f:
        chars = f.readline().strip().split()

        for line in f:
            ch1, *scores = line.strip().split()

            for i, score in enumerate(scores):
                scoring_matrix[(ch1, chars[i])] = int(score)

    return scoring_matrix


def read_seq_from_file(seq_file):
    seq_id = ''
    seq = ''
    with open(seq_file) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
            else:
                seq += line.strip()
    return seq_id, seq


def read_seqs(seqs):
    dict_seqs = {}
    for seq in seqs:
        seq_id, seq = read_seq_from_file(seq)
        dict_seqs[seq_id] = seq
    return dict_seqs


def get_from_folder_and_insert():
    seqs = {}

    id1, seq1 = read_seq_from_file('./genomes/DQ182595.1_SARS.fasta')
    seqs[id1] = seq1

    id2, seq2 = read_seq_from_file('./genomes/DQ648857.1_BAT.fasta')
    seqs[id2] = seq2

    id3, seq3 = read_seq_from_file('./genomes/EPI_ISL_402131_BAT.fasta')
    seqs[id3] = seq3

    id4, seq4 = read_seq_from_file('./genomes/EPI_ISL_404253_SARS-CoV-2.fasta')
    seqs[id4] = seq4

    id5, seq5 = read_seq_from_file('./genomes/JX869059.2_MERS.fasta')
    seqs[id5] = seq5

    id6, seq6 = read_seq_from_file('./genomes/JX993987.1_BAT.fasta')
    seqs[id6] = seq6

    id7, seq7 = read_seq_from_file('./genomes/KT368829.1_MERS.fasta')
    seqs[id7] = seq7

    id8, seq8 = read_seq_from_file('./genomes/NC_004718.3_SARS.fasta')
    seqs[id8] = seq8

    id9, seq9 = read_seq_from_file('./genomes/NC_045512.2_SARS-CoV-2.fasta')
    seqs[id9] = seq9
    return seqs


def build_db(db):
    db_dict = {}

    for i in range(0, len(db) - k + 1):
        if db[i:i + k] in db_dict:
            db_dict[db[i:i + k]].append(i)
        else:
            db_dict[db[i:i + k]] = [i]

    return db_dict


def align(seq1, seq2, scoring_matrix):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be with the same length")
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be with the same length")
    score = 0
    for i in range(0, len(seq1)):
        score = score + scoring_matrix[seq1[i], seq2[i]]

    return score


def find_neighbors(kmer, scoring_matrix, alphabet):
    neighbors = []
    max_score = align(kmer, kmer, scoring_matrix)

    if max_score >= T:
        find_neighbors_rec(kmer, kmer, 0, max_score, alphabet, neighbors, scoring_matrix)

    return neighbors


def find_neighbors_rec(kmer, neighbor, pos, curr_score, alphabet, neighbors, scoring_matrix):
    if len(kmer) == pos:
        neighbors.append(neighbor)
    else:
        for i in range(0, len(alphabet)):
            neighbor = neighbor[:pos] + alphabet[i] + neighbor[pos + 1:]
            score = curr_score - scoring_matrix[kmer[pos], kmer[pos]] + scoring_matrix[kmer[pos], alphabet[i]]
            if score >= T:
                find_neighbors_rec(kmer, neighbor, pos + 1, score, alphabet, neighbors, scoring_matrix)


def get_hsps(query, db_dict, scoring_matrix, alphabet):
    hsps = []

    for i in range(0, len(query) - k + 1):
        kmer = query[i:i + k]
        neighbors = find_neighbors(kmer, scoring_matrix, alphabet)
        for val in neighbors:
            if val in db_dict:
                score = align(kmer, val, scoring_matrix)
                for j in db_dict[val]:
                    hsps.append(HSP(seq2_start=i, seq2_end=i + k - 1, seq1_start=j, seq1_end=j + k - 1, score=score))

    return hsps


def find_HSP(seq1, seq2, matrix, alphabet):
    db_dict = build_db(seq2)
    return get_hsps(seq1, db_dict, matrix, alphabet)


def extend_left(query, db, hsp, scoring_matrix):
    """returns the left extension with the maximal score"""

    msp = copy(hsp)

    maxScore = msp.score
    mspMark = msp.seq1_start
    dbMark = msp.seq2_start
    good = 0

    sum_so_far = msp.score

    while (mspMark > 0) and (dbMark > 0) and (good == 0):
        mspMark = mspMark - 1
        dbMark = dbMark - 1

        sum_so_far = sum_so_far + scoring_matrix[query[mspMark], db[dbMark]]

        scoreCheck = maxScore - X

        if sum_so_far >= maxScore:
            maxScore = sum_so_far
            msp.score = maxScore
            msp.seq1_start = mspMark
            msp.seq2_start = dbMark
        elif sum_so_far < scoreCheck:
            good = 1

    return msp


def extend_right(query, db, hsp, scoring_matrix):
    """returns the right extension with the maximal score"""

    msp = copy(hsp)

    maxScore = msp.score
    mspMark = msp.seq1_end
    dbMark = msp.seq2_end
    good = 0
    sum_so_far = msp.score

    while (mspMark < len(query) - 1) and (dbMark < len(db) - 1) and (good == 0):
        mspMark = mspMark + 1
        dbMark = dbMark + 1

        sum_so_far = sum_so_far + scoring_matrix[query[mspMark], db[dbMark]]

        scoreCheck = maxScore - X

        if sum_so_far > maxScore:

            maxScore = sum_so_far
            msp.score = maxScore
            msp.seq1_end = mspMark
            msp.seq2_end = dbMark
        elif sum_so_far < scoreCheck:
            good = 1
    return msp


def extend_hsp(query, db, hsp, scoring_matrix):
    msp_left = extend_left(query, db, hsp, scoring_matrix)
    msp = extend_right(query, db, msp_left, scoring_matrix)

    return msp


def extend_HSP_to_MSP(hsps, seq1, seq2, matrix, R):  # order in adding? keeping the last.
    msp_by_diagonals = {}
    msps = []

    for hsp in hsps:
        if hsp.diagonal() not in msp_by_diagonals.keys() \
                or (hsp.seq1_start > msp_by_diagonals[hsp.diagonal()].seq1_end and abs(hsp.diagonal()) < R):
            msp = extend_hsp(seq2, seq1, hsp, matrix)  # error out of bound
            msp_by_diagonals[hsp.diagonal()] = msp
            msps.append(msp)
    return msps


def need_to_add(seq1_end, seq2_end, seq1_start, seq2_start):
    return seq1_end < seq1_start and seq2_end < seq2_start


def update_graphs(graph, msps):
    for msp1 in msps:
        for msp2 in msps:
            if not msp1.__eq__(msp2):
                if need_to_add(msp1.seq1_end, msp1.seq2_end, msp2.seq1_start, msp2.seq2_start):
                    graph.add_edge(msp1, msp2, weight=msp1.score)


def add_end_node(graph):
    end_node = HSP(0, 0, 0, 0, -1)  # maybe need to change
    graph.add_node(end_node)
    for node in graph.nodes:
        if not node.__eq__(end_node):
            graph.add_edge(node, end_node, weight=node.score)


def add_nodes(graph, msps):
    for msp in msps:
        graph.add_node(msp)


def create_graphs(msps):  # check!

    graph = nx.DiGraph()
    add_nodes(graph, msps)  # adding all msps
    update_graphs(graph, msps)  # adding all edges
    add_end_node(graph)  # adding the last node and edges
    return graph


def find_heaviest_paths(graph):
    path = nx.dag_longest_path(graph, weight='weight')
    score = 0
    for hsp in path:
        score += hsp.score
    return score


# used to print fot the notebook
def count_msp(msps_by_pairs):
    for pair, msps in msps_by_pairs.items():
        print(pair, len(msps))
        print("----------------------")


def update_output_file(seq1_name, seq2_name, score):
    output_file = open("scores.txt", "a+")
    output_line = f"{seq1_name}\t{seq2_name}\t{score}\t\n"
    output_file.write(output_line)
    output_file.close()


def parse_matrix(path):
    last_name = ''

    seqs_names = []
    to_add_line = []
    big_matrix = []
    not_first = False
    need_to_add_name = True

    with open(path) as f:
        for line in f:
            splited = line.strip('\n').split('\t')
            name = splited[0]
            if not name == last_name:
                last_name = name
                if need_to_add_name and not seqs_names.__contains__(splited[1]):
                    to_add_line.append(0)  # adding the compare to itself
                    seqs_names.append(name)
                if not_first:
                    to_add_line.reverse()
                    big_matrix.append(to_add_line)
                    to_add_line = []
                    to_add_line.append(0)  # adding the compare to itself
                    need_to_add_name = False
                not_first = True
            if need_to_add_name and not seqs_names.__contains__(splited[1]):
                seqs_names.append(splited[1])
            to_add_line.append(int(splited[2]))

        to_add_line.reverse()
        big_matrix.append(to_add_line)
        big_matrix.append([0])

        big_matrix.reverse()
        seqs_names.reverse()  # chek if correct

        print(big_matrix)
        # print (np.array(big_matrix))
        return big_matrix, seqs_names


# calculate R - if we consider the length of the sequences, we can assume the best match won't be further than
# 85% from the average of the sequences.

def get_R(seqs_dict):
    len_of_all_seqs = 0
    num_seqs = 0
    for seq_id, seq in seqs_dict.items():
        num_seqs += 1
        len_of_all_seqs += len(seq)
    return (len_of_all_seqs / num_seqs) * 0.85


def main():
    alphabet = get_alphabet(sys.argv[1])
    matrix = read_scoring_matrix(sys.argv[1])
    seqs = sys.argv[2:]
    seqs_value = read_seqs(seqs)

    R = get_R(seqs_value)

    i = j = 0

    # compute for all the pairs of seqs
    for id1, seq1 in seqs_value.items():
        for id2, seq2 in seqs_value.items():
            if i < j:
                start_time = datetime.now()
                # print(i, j)

                hsps = find_HSP(seq1, seq2, matrix, alphabet)  # section 1 - done
                # print("after HSP")
                msps = extend_HSP_to_MSP(hsps, seq1, seq2, matrix, R)  # section 2
                # print("after MSP")

                seqs_graph = create_graphs(msps)  # section 3

                heaviest_path_score = find_heaviest_paths(seqs_graph)  # section 4

                delta_time = datetime.now() - start_time
                update_output_file(id1, id2, heaviest_path_score)
                print(delta_time)
            j += 1
        j = 0
        i += 1


if __name__ == '__main__':
    main()
