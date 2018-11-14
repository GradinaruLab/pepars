from __future__ import print_function, division
import argparse
import sys
import re

def mod_distance(s1, s2):
    """Given two same-length strings, return the number of characters that are different between s1 and s2."""
    return sum(1 if a != b else 0 for a, b in zip(s1, s2))


def dist_is_1(s1, s2):
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
            if count == 2:
                return False
    return count != 0


def try_int(x):
    if x.strip() == "":
        return 0

    try:
        return int(x.strip())
    except:
        return x


def read_csv(fn):
    """Return a list of (string, count) tuples from the csv file."""
    with open(fn) as f:
        entries = []
        for line in f:
            toks = [x.strip() for x in line.split(",")]
            if len(toks) <= 1:
                continue
            entries.append([try_int(x) for x in toks])
    return entries

ALTERNATIVES = {"A": frozenset("CGT"),
                "C": frozenset("AGT"),
                "G": frozenset("ACT"),
                "T": frozenset("ACG")}


def string_neighbors(s):
    for i in range(len(s)):
        for letter in ALTERNATIVES[s[i]]:
            yield s[:i] + letter + s[i+1:]


def all_neighbors(entries, sc):
    neighbors = {}
    for i, e in reversed(list(enumerate(entries))):
        for s in string_neighbors(e[sc]):
            neighbors[s] = i
    return neighbors


def filter_similar(entries, sc, cc, nc, top_n):
    """Given a set of entries, return a filtered version where entries in
bottom_pct differing by one-letter from top_pct are removed."""
    bottom_n = len(entries) - top_n

    print("Filtering bottom {} rows based on top {} rows of column {} (out of {} rows total).".format(bottom_n, top_n, cc, len(entries)), file=sys.stderr)

    new_entries = entries[:-bottom_n]
    neighbors = all_neighbors(entries[:top_n], sc)

    removed_entries = []
    for e in entries[-bottom_n:]:
        if e[sc] in neighbors:
            n = neighbors[e[sc]]
            parent = entries[n]
            removed_entries.append((e, n, parent[nc], parent[sc], cc))
        else:
            new_entries.append(e)

    return new_entries, removed_entries


def parse_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--top-n", "--top", action="store", default=None,
                        help="Filter based on the top number of rows")

    parser.add_argument("--count-columns", "-c", action="store", default="3",
                        help="Columns containing the counts to sort by")
    parser.add_argument("--name-column", "-n", action="store", default=2, type=int,
                        help="Column storing the name to reference")
    parser.add_argument("--string-column", "-f", action="store", default=1, type=int,
                        help="Column containing the strings to differ by")

    parser.add_argument("input_csv", action="store",
                        help="input file to read")

    return parser.parse_args()


def main():
    args = parse_args()
    entries = read_csv(args.input_csv)
    count_columns = [int(x) for x in re.split(r"[, ]+", args.count_columns)]
    tops = [int(x) for x in re.split(r"[, ]+", args.top_n)]
    if len(tops) != len(count_columns):
        if len(tops) == 1:
            tops = tops * len(count_columns)
        else:
            print("ERROR: Number of --top specified doesn't match number of columns", file=sys.stderr)
            return 1

    removed = []
    nc = args.name_column - 1
    sc = args.string_column - 1
    for cc, t in zip(count_columns, tops):
        entries.sort(key=lambda x: x[cc-1], reverse=True)

        filtered, new_removed = filter_similar(entries, sc, cc, nc, t)
        print("Removed {} entries based on Column {}".format(len(new_removed), cc), file=sys.stderr)
        entries = filtered

        removed.extend(new_removed)

    for e in entries:
        print(",".join(str(x) for x in e))
    print("*")

    print("{} Removed Entries".format(len(removed)))
    for e, parent_row, orig_name, orig_string, col in removed:
        print(",".join(str(x) for x in e) +
              ",removed due to column {} and similar to {} (Row {} - {})".format(col, orig_string, parent_row + 1, orig_name))

if __name__ == "__main__":
    main()
