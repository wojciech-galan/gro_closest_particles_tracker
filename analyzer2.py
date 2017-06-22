#! /usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import math
import argparse
import itertools
import numpy as np

frame_expression = re.compile('^REMARK.+ENDMDL$', re.DOTALL|re.MULTILINE)
model_expression = re.compile('^MODEL', re.MULTILINE)
ter_expression = re.compile('^TER', re.MULTILINE)


class AtomProperties(object):

    def __init__(self, num, atom, particle, coordinates):
        super(AtomProperties, self).__init__()
        self.number = num
        self.atom = atom
        self.particle = particle
        self.coordinates = coordinates

    def __repr__(self):
        return 'AtomProperties({number}, {atom}, {particle}, {coordinates!r})'.format(**self.__dict__)


def get_frames(f_name, bytes):
    '''Zwraca kolejne ramki'''
    with open(f_name) as f:
        # http://stackoverflow.com/questions/32040009/find-out-if-file-pointer-is-at-eof-in-python
        f.seek(0, os.SEEK_END)
        n = f.tell()
        f.seek(0)
        line = ''
        x = None
        while f.tell()<n:
            text = line
            while not x and f.tell()<n:
                line = f.read(200000)
                text += line
                x = re.findall(frame_expression, text)
            assert len(x)==1
            yield x[0].strip()
            x = None

def parse_frame(frame_text, particle_):
    '''Zwraca współrzędne atomów w analizowanej cząsteczce, oraz atomów tlenu z wody'''
    particle_atoms, water_oxygens = [], []
    model_part = re.search(model_expression, frame_text)
    _, model, rest = frame_text[model_part.start():].split(None, 2)
    ter_part = re.search(ter_expression, rest)
    new_rest = rest[:ter_part.start()].strip()
    # teraz każda z tych linii powinna zaczynac sie od ATOM
    for line in new_rest.split(os.linesep):
        _, num, atom, particle, _, x, y, z, _, _, _ = line.split()
        if particle == particle_:
            particle_atoms.append(AtomProperties(num, atom, particle, (float(x), float(y), float(z))))
        elif particle == 'SOL' and atom == 'OW':
            water_oxygens.append(AtomProperties(num, atom, particle, (float(x), float(y), float(z))))
    return particle_atoms, water_oxygens

def distance((x1, y1, z1), (x2, y2, z2)):
    '''liczy odległość euklidesową'''
    return math.sqrt((x1-x2) ** 2 + (y1-y2) ** 2 + (z1-z2) ** 2)

def find_closest(cholesterol_atoms, water_oxygens, dist):
    """
    Zwraca numery atomów spośród water_oxygens, które sa odległe o co najwyżej dist od któregokolwiek
    z cholesterol_atoms
    :return: set
    """
    out_nums = set()
    for oxygen in water_oxygens:
        for atom in cholesterol_atoms:
            if distance(oxygen.coordinates, atom.coordinates) <= dist:
                out_nums.add(oxygen.number)
                break
    return out_nums

def how_many_left(cholesterol_atoms, water_oxygens, atom_numbers, dist):
    """
    Sprawdza, jak wiele z cząsteczek wody (konkretnie atomów tlenu) o numerach sposród atom_numbers pozostało ciagle
    w oległości nie większej niż dist od cholesterolu
    :return: int
    """
    num = 0
    for oxygen in water_oxygens:
        if oxygen.number in atom_numbers:
            for atom in cholesterol_atoms:
                if distance(oxygen.coordinates, atom.coordinates) <= dist:
                    num += 1
                    break
    return num

def analysis(infile, frame_count, particle, distance, separate):
    i = 0
    how_many = []
    if separate:
        how_many_separate = []
    for frame in get_frames(infile, 200000):
        particle_atoms, water_oxygens = parse_frame(frame, particle)
        if separate:
            separate_atoms = []
            for atom in particle_atoms:
                if atom.atom in separate:
                    separate_atoms.append(atom)
                    particle_atoms.remove(atom)
        if not i % frame_count:
            atom_nums = find_closest(particle_atoms, water_oxygens, distance)
            how_many.append(len(atom_nums))
            if separate:
                separate_atom_nums = find_closest(separate_atoms, water_oxygens, distance)
                how_many_separate.append(len(separate_atom_nums))
        else:
            how_many.append(how_many_left(particle_atoms, water_oxygens, atom_nums, distance))
            if separate:
                how_many_separate.append(how_many_left(separate_atoms, water_oxygens, separate_atom_nums, distance))
        i += 1
    ar = np.array([how_many[x:x + frame_count] for x in range(0, len(how_many) / frame_count * frame_count, frame_count)])
    if separate:
        sep_ar = np.array([how_many_separate[x:x + frame_count] for x in range(0, len(how_many_separate) /
                                                                               frame_count * frame_count, frame_count)])
        return np.array([np.mean(ar, axis=0), np.std(ar, axis=0), np.mean(sep_ar, axis=0), np.std(sep_ar, axis=0)])
    return np.array([np.mean(ar, axis=0), np.std(ar, axis=0)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("frames", help='plik z ramkami')
    parser.add_argument('outfile', help='plik wynikowy')
    parser.add_argument('frame_count', help='ile ramek. Np. 100 oznacza "ramkę 0" + 99 kolejnych, w których \
    obserwujesz zachowanie wody', type=int)
    parser.add_argument('--distance', '-d', help='odległość od cholesterolu w Angstremach', type=int, default=4)
    parser.add_argument('--particle', '-p', help='cząsteczka, którą chcesz analizowac', type=str, default='CHL')
    parser.add_argument('--separate', '-s',
                        help='''atomy, które chcesz analizować osobno od reszty analizowanej cząsteczki.
                        Oddzielone przecinkami, np: HO3,O3''',
                        type=str, default='')
    args = parser.parse_args()
    res = analysis(args.frames, args.frame_count, args.particle, args.distance, tuple(args.separate.split(',')))
    print res
    np.savetxt(args.outfile, res.transpose(), delimiter=',')

    # i = 0
    # dzielnik = 100#25
    # how_many = []
    # for frame in get_frames('../data/ramki.pdb', 200000):
    #     cholesterol_atoms, water_oxygens =  parse_frame(frame)
    #     if not i%dzielnik:
    #         atom_nums = find_closest(cholesterol_atoms, water_oxygens, 4)
    #         how_many.append(len(atom_nums))
    #         #print [how_many[x:x+dzielnik] for x in range(0, len(how_many), dzielnik)]
    #     else:
    #         how_many.append(how_many_left(cholesterol_atoms, water_oxygens, atom_nums, 4))
    #     i += 1
    #
    # print len(how_many)/float(dzielnik)
    # ar = np.array([how_many[x:x + dzielnik] for x in range(0, len(how_many)/dzielnik*dzielnik, dzielnik)])
    # print np.mean(ar, axis=0)




