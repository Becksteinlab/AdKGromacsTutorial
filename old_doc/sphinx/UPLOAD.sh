#!/bin/bash
PRACTICAL=13
HTML=AdKTutorial

make html
(cd .. && tar -jcvf AdKTutorial.tar.bz2 AdKTutorial/ )

rsync -avP --delete _build/html/ becksteinlab:public_html/pages/courses/2013/SimBioNano/${PRACTICAL}/${HTML}
cput $PRACTICAL ../AdKTutorial.tar.bz2 