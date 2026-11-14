"""Create SHA256SUMS for one release's wheel and source archive."""
import argparse
import hashlib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    wheels = list(args.directory.glob('*.whl'))
    sources = list(args.directory.glob('*.tar.gz'))
    if len(wheels) != 1 or len(sources) != 1:
        parser.error('Expected exactly one wheel and one source archive in directory.')
    lines = []
    for path in sorted(wheels + sources):
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        lines.append(f'{digest.hexdigest()}  {path.name}\n')
    (args.directory / 'SHA256SUMS').write_text(''.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
