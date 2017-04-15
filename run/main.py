#!/usr/bin/env python3
# -*- coding: utf8 -*-


import argparse
import noahmp_config


DEFAULT_NAMELIST_FILE = 'case.nml'


def main(nmlfile):
    cfg = noahmp_config.Config(nmlfile)
    print(dir(cfg))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Noah-MP Land Surface Model')
    parser.add_argument('nmlfile', nargs='?', type=str,
                        default=DEFAULT_NAMELIST_FILE,
                        help='configuration file')
    args = parser.parse_args()
    main(args.nmlfile)
