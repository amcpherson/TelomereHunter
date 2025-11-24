from setuptools import setup, find_packages

setup(
    name="telomerehunter",
    version="1.1.0",
    description="Estimation of Telomere Content from WGS Data",
    long_description=open("README.md").read() if True else "",
    long_description_content_type="text/x-rst",
    author="Lars Feuerbach, Philip Ginsbach, Lina Sieverling",
    author_email="l.sieverling@dkfz-heidelberg.de",
    url="https://www.dkfz.de/en/applied-bioinformatics/telomerehunter/telomerehunter.html",
    license="GPL",
    packages=find_packages(),
    include_package_data=True,
    keywords=["telomere", "content", "read", "NGS", "WGS", "tumor", "control"],
    scripts=[
        "scripts/telomerehunter",
    ],
    install_requires=[
        "PyPDF2",
        "numpy",
        "pysam",
    ],
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
)

