def remove_from_file(chimeras, file):
    with open(file, "r+") as f:
        d = f.readlines()
        f.seek(0)
        count = 0
        while count < len(d):
            if any(chimera in d[count] for chimera in chimeras):
                count += 4
                continue
            else:
                f.write(d[count])
            count += 1
        f.truncate()