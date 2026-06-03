#!/usr/bin/env python3
"""
rebind_pdf.py  —  Bind all JPGs in the JPGs/ directory into a single PDF.

Usage
-----
  python rebind_pdf.py                      # default: JPGs/ → JPGs/all_examples.pdf
  python rebind_pdf.py --dir MyFigs         # different input directory
  python rebind_pdf.py --out output.pdf     # different output filename
"""

import argparse
import glob
import os
import sys


def bind(jpg_dir, output_pdf):
    jpgs = sorted(glob.glob(os.path.join(jpg_dir, '*.jpg')))

    if not jpgs:
        sys.exit(f"No JPG files found in '{jpg_dir}'.")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        sys.exit("matplotlib is required.  Install with:  pip install matplotlib")

    with PdfPages(output_pdf) as pdf:
        for jpg in jpgs:
            img = mpimg.imread(jpg)
            h, w = img.shape[:2]
            fig, ax = plt.subplots(figsize=(w / 100, h / 100))
            ax.imshow(img)
            ax.axis('off')
            plt.tight_layout(pad=0)
            pdf.savefig(fig, bbox_inches='tight', dpi=150)
            plt.close(fig)
            print(f"  {os.path.basename(jpg)}")

    print(f"\nPDF saved: {output_pdf}  ({len(jpgs)} pages)")


def main():
    parser = argparse.ArgumentParser(
        description="Bind all JPGs in a directory into a single PDF.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--dir", default="JPGs",
                        help="directory containing the JPG files (default: JPGs)")
    parser.add_argument("--out", default=None,
                        help="output PDF path (default: <dir>/all_examples.pdf)")
    args = parser.parse_args()

    output = args.out or os.path.join(args.dir, "all_examples.pdf")
    bind(args.dir, output)


if __name__ == "__main__":
    main()
