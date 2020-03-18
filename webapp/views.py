# The MIT License (MIT)
# Copyright (c) 2019 Thomas Huetter, Manuel Kocher.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# import needed django modules
from django.shortcuts import render
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponseRedirect
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponse
from wsgiref.util import FileWrapper
import json

# import designate.py
from webapp.designate import *

# import modules to compute alignments and trees
from Bio import Phylo
from Bio import AlignIO

# import temporary files and zip files
import tempfile
from zipfile import ZipFile
import os
from django.conf import settings

# imports for testing and debug
from io import StringIO
import time
import sys
import logging

# logging config for debugging only
logging.basicConfig(format='%(asctime)s %(message)s')
logger=logging.getLogger()
logger.setLevel(logging.DEBUG)

'''
    @param request: httprequest of webpage home.html
    @return HttpResponseRedirect: redirect to url /files
        -   when the button "get started" is clicked, the user will be
            redirected to /files.
    @return render: render the webpage home.html
        -   the webpage fileInput.html renders if the next button is not clicked
'''
def home(request):

    # flush possibly still existing older session data
    request.session.flush()

    if request.method == 'POST':
        return HttpResponseRedirect('/files')

    return render(request, 'home.html')

'''
    @param request: httprequest of webpage home.html
    @return HttpResponseRedirect: redirect to url /tree
        -   when the button "next" is clicked, the alignments and the newick
            file is computed. if it is successful, the alignment and tree are
            saved in the session and the user is redirected to /tree.
    @return render: render the webpage home.html
        -   the webpage fileInput.html renders if the next button is not clicked
        -   if the files cannot be computed, then an error message is displayed
            and the user stays at fileInput.html
'''
def file_input(request):

    # flush possibly still existing older session data
    request.session.flush()

    # compute data if button is clicked (e.g. if method is POST)
    if request.method == 'POST':
        try:
            if request.POST.get('exampleFile') == 'on':
                # use standard example files (checkbox is clicked)
                tree_file = open(settings.BASE_DIR + "/example_input/example.fas.treefile", "rb").read()
                align_file = open(settings.BASE_DIR + "/example_input/example.fas", "rb").read()
            else:
                # get necessary alignment and tree files
                if "newick" in request.FILES:
                    tree_file = request.FILES['newick'].read()
                else:
                    tree_file = False
                align_file = request.FILES['fasta'].read()

            # get k value from input
            request.session['k_value'] = request.POST.get('kValue')

            if request.POST.get('considerGaps') == 'on':
                request.session['consider_gaps'] = True
            else:
                request.session['consider_gaps'] = False

            # write to temporary files and use AlignIO and Phylo
            with tempfile.NamedTemporaryFile() as temp:
                temp.write(align_file)
                temp.seek(0)
                request.session['alignments'] = AlignIO.read(temp.name, "fasta")
            # delete tempfile
            temp.close()

            if tree_file == False:
                return HttpResponseRedirect('/names')

            else:
                with tempfile.NamedTemporaryFile() as temp:
                    temp.write(tree_file)
                    temp.seek(0)
                    request.session['phylo_tree'] = Phylo.read(temp.name, "newick")
                # delete tempfile
                temp.close()

            # read newick file as string
            request.session['newick_tree'] = tree_file.splitlines()

            # assign preorder id
            request.session['phylo_tree'] = assign_pre_id_to_inner_nodes(request.session.get('phylo_tree'))

            # TODO: - function to check if uploaded files are really fasta and
            #         newick files
            #       - rewrite file input to match thomas program.py functions
            #         (files are already computed using alignIO and phylo)

            return HttpResponseRedirect('/tree')

        # if an error occures, show the error message and flush the session data
        except Exception as e:
            error_message = str(e) + ". Please return to the main page."
            request.session.flush()

            return render(request, 'fileInput.html', {'error_message':error_message})

    return render(request, 'fileInput.html')

'''
    @param request: httprequest of webpage tree.html
    @return HttpResponseRedirect: redirect to url /settings
        -   if the next button is clicked and enough nodes were selected,
            the user is redirected to /settings
    @return render: render the webpage tree.html
        -   the webpage tree.html renders if the next button is not clicked
        -   if an error occures, an specific error message will appear
'''
def tree(request):

    # compute data if button is clicked (e.g. if method is POST)
    if request.method == 'POST':
        try:

            # get ids of selected nodes
            id = request.POST.get('hiddenElement')

            # raise error if last or first node is empty
            if id[0] == "," or id[-1] == ",":
                raise ValueError()

            # go through characters and convert nodes into a list
            lst = []
            tmp = ""
            for c in id:
                if c == ',':
                    lst.append(tmp)
                    tmp = ""
                    continue
                tmp = tmp + c
            lst.append(tmp)

            # if less or more than 2 nodes are selected, raise an error
            if not len(lst) == 2:
                raise ValueError()

            # save node ids in session
            request.session['id0'] = lst[0]
            request.session['id1'] = lst[1]

            return HttpResponseRedirect('/results')

        # throw error if not enough nodes are selected
        # show tree again for convenience
        except ValueError as e:
            error_message = "Too few nodes were selected. Please select 2 nodes (Query and Reference group)"
            num_leafs = len(get_all_terminals_of_inner_node(request.session.get('phylo_tree'), '0'))
            depth = get_depth(request.session.get('phylo_tree'))

            return render(request, 'tree.html', {"error_message": error_message, 'newickTree': request.session.get('newick_tree'), 'numLeafs': num_leafs, 'depth': depth})

        # standard error message, flush session memory
        except Exception as e:
            error_message = str(e) + ". Please return to the main page."
            request.session.flush()

            return render(request, 'tree.html', {"error_message": error_message})

    try:
        # get number of leafs for width calc and depth for height calc
        num_leafs = len(get_all_terminals_of_inner_node(request.session.get('phylo_tree'), '0'))
        depth = get_depth(request.session.get('phylo_tree'))

        return render(request, 'tree.html', {'newickTree': request.session.get('newick_tree'), 'numLeafs': num_leafs, 'depth': depth})

    # this error message appears when the user jumps directly to tree without
    # uploading any file whatsoever
    except AttributeError as e:
        error_message = str(e) + ". Did you upload files in the right format? Please try again."
        request.session.flush()

        return render(request, 'tree.html', {"error_message": error_message})

    # standard error message
    except Exception as e:
        error_message = str(e) + ". Please return to the main page."
        request.session.flush()

        return render(request, 'tree.html', {"error_message": error_message})

def names(request):

    if request.method == 'POST':
        # get ids of selected nodes
        request.session['qGroup'] = [x.strip() for x in request.POST.get('qGroup').split(',')]
        request.session['rGroup'] = [x.strip() for x in request.POST.get('rGroup').split(',')]

        request.session['id0'] = request.session['rGroup']
        request.session['id1'] = request.session['qGroup']


        request.session['phylo_tree'] = None


        return HttpResponseRedirect('/results')

    try:
        return render(request, 'tree2.html', {'alignments': [x.name for x in request.session['alignments']]})

    except Exception as e:
        error_message = "An error occured. Did you upload a .fas file in the right format? Please return to the main page."
        request.session.flush()

        return render(request, 'tree2.html', {"error_message": error_message})

'''
    @param request: httprequest of webpage results.html
    @return render: render the webpage settings.html
        -   the webpage results.html renders with calculated results
        -   if an error occures, and expection is thrown and
            an specific error message will appear, python errormessage otherwise
'''
def results(request):
    try:
        if request.session['phylo_tree'] != None:
            # compute alignments of both groups
            alignments_ref_group = get_alignments_of_group(
                get_all_terminals_of_inner_node(request.session.get('phylo_tree'), request.session.get('id0')), request.session.get('alignments'))
            alignments_query_group = get_alignments_of_group(
                get_all_terminals_of_inner_node(request.session.get('phylo_tree'), request.session.get('id1')), request.session.get('alignments'))

        else:
            alignments_ref_group = get_groups_from_list(request.session['rGroup'], request.session.get('alignments'))
            alignments_query_group = get_groups_from_list(request.session['qGroup'], request.session.get('alignments'))


        # remove possible duplicates in second group
        alignments_ref_group = list_difference(alignments_query_group,
                                                    alignments_ref_group)

        # compute d-power, q-rank and r-rank with individual nucleotide and combined nucleotide analysis
        ranking = signature_character_detection(request.session.get('alignments'), request.session.get('phylo_tree'), request.session.get('id1'),
                                                request.session.get('id0'), int(request.session.get('k_value')), None,False,
                                                False, request.session.get('consider_gaps'))

        request.session['ranking'] = ranking


        # Missing group id in function call computes the shannon entropy for the full alignment.
        entropy_list_temp = shannon_entropy_analysis(request.session.get('alignments'), request.session.get('phylo_tree'))


        request.session['entropy_list_temp'] = entropy_list_temp

        # extract the exact and average entropy from the analysis (for chart.js)
        entropy_list = []
        average_list = []
        for element in entropy_list_temp:
            entropy_list.append(element[1])
            average_list.append(element[2])

        # since there is no easy way to handle multidimensional arrays in
        # the template, a list with binary, asymmetric, noisy must be constructed
        # where [0] is binary, [1] is asymmetric and [2] is noisy
        diversity_list = [[], [], []]
        nucleotide_list = [[], []]
        id = 0
        combined_asymmetric = 0
        while id < len(ranking):
            if len(ranking[id]) > 4:
                if ranking[id][4] == 'binary':
                    diversity_list[0].append(ranking[id][0])

                    # nucleotide_list[0] contains binary positions with nucleotide letter
                    nucleotide_list[0].append([ranking[id][0], ranking[id][2][1][0].upper()])
                if ranking[id][4] == 'asymmetric':
                    if isinstance(ranking[id][0], int):
                        diversity_list[1].append(ranking[id][0])

                        # nucleotide_list[1] contains asymmetric positions with nucleotide letter
                        nucleotide_list[1].append([ranking[id][0], ranking[id][2][1][0].upper()])
                    else:
                        combined_asymmetric = combined_asymmetric + 1
                if ranking[id][4] == 'noisy':
                    diversity_list[2].append(ranking[id][0])
            id = id + 1


        # sort the lists by position ASC to save searches later on
        diversity_list[0] = sorted(diversity_list[0])
        diversity_list[1] = sorted(diversity_list[1])
        diversity_list[2] = sorted(diversity_list[2])

        nucleotide_list[0] = sorted(nucleotide_list[0], key=lambda x:x[0])
        nucleotide_list[1] = sorted(nucleotide_list[1], key=lambda x:x[0])

        # convert python into json object for compatibility with template
        for x in nucleotide_list[0]:
            x = [x[0], json.dumps(x[1])]
        for x in nucleotide_list[1]:
            x = [x[0], json.dumps(x[1])]

        # create a string for the files with group ids and session key for uniqueness
        ids_key = str(str(request.session.get('id0')) + "_" + str(request.session.get('id1')) + "_" + request.session.session_key)

        request.session['ids_key'] = ids_key

        return render(request, 'results.html', {'alignment_reference': alignments_ref_group,
                'alignment_query': alignments_query_group, 'entropy_list':entropy_list,
                'average_list': average_list,
                'diversity_list': diversity_list, 'ranking': ranking,
                'nucleotide_list': nucleotide_list, 'ids': ids_key,
                'combined_asymmetric': combined_asymmetric, 'k_window': request.session.get('k_value')
                })

    # throw error if no values to compute are available
    except AttributeError as e:
        error_message = "No input data is available. Please return to the main page."
        request.session.flush()

        return render(request, 'results.html', {'error_message': error_message})

    # throw error if selected values could not be found
    except SystemError as e:
        error_message = "The selected names could not be found in the input data."
        request.session.flush()

        return render(request, 'results.html', {"error_message": error_message})

    # standard error message
    except Exception as e:
        error_message = str(e) + ". Please return to the main page."
        request.session.flush()

        return render(request, 'results.html', {"error_message": error_message})

'''
    @param request: httprequest of webpage download_zip.html
    @return response: return a HttpResponse
        -   return a response which asks the user if he wants to download the .zip file
'''
def download_zip(request):

    # prepare the files for download
    create_zip(request, request.session.get('ranking'), request.session.get('entropy_list_temp'), request.session.get('ids_key'))

    # TODO: Is it possible to create the zip only in this case? Do we have the ranking and the entropy here?
    ids = str(request.session.get('id0')) + "_" + str(request.session.get('id1'))
    filename = settings.MEDIA_ROOT + "DeSigNate_Groups_" + time.asctime().replace(' ', '_').replace(':', '_') + "_" + request.session.session_key + ".zip"

    response = HttpResponse(open(filename, 'rb'), content_type='application/zip')
    response['Content-Length'] = os.path.getsize(filename)
    response['Content-Disposition'] = 'attachment; filename=%s' % str('DeSigNate_Results_' + time.asctime().replace(' ', '_').replace(':', '_') + '.zip')
    return response

def contact(request):
    return render(request, 'contact.html')

def privacy(request):
    return render(request, 'privacy.html')

def impressum(request):
    return render(request, 'impressum.html')

'''
    @param ranking: list with position, d-power, q-rank and r-rank
    @param entropy: list of entropy values create by shannon_entropy_analysis()
    @param ids:     string consisting of the query group id, the reference group id, and the session_key
'''
def create_zip(request, ranking, entropy, ids):
    ranking_filepath = settings.MEDIA_ROOT + "ranking.csv"
    entropy_filepath = settings.MEDIA_ROOT + "entropy.csv"
    zip_filepath = settings.MEDIA_ROOT + "DeSigNate_Groups_" + time.asctime().replace(' ', '_').replace(':', '_') + '_' + request.session.session_key + ".zip"

    write_csv_file(ranking_filepath, ["position", "discriminative_power",
        "query_rank", "reference_rank"], ranking)
    write_csv_file(entropy_filepath, ["position", "shannon_entropy", "average"], entropy)

    files_to_zip = [entropy_filepath, ranking_filepath]

    with ZipFile(zip_filepath, 'w') as zip:
        # writing each file one by one
        for i in range(len(files_to_zip)):
            zip.write(files_to_zip[i], os.path.basename(files_to_zip[i]))

    os.chmod(zip_filepath, 0o777)
    os.remove(ranking_filepath)
    os.remove(entropy_filepath)
