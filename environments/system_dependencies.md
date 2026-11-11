# System-level dependencies

The following tools were installed system-wide on Ubuntu
and were not part of the exported micromamba environments.

| Tool | Version | Executable |
|---|---|---|
| SeqKit | 2.3.0 | /usr/bin/seqkit |
| Samtools | 1.19.2 | /usr/bin/samtools |

Installation on Ubuntu:

sudo apt update
sudo apt install seqkit samtools

Verify installations:

seqkit version
samtools --version
