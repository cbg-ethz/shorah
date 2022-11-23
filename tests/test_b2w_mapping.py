from array import array
import pytest
from cigar import Cigar
from shorah import b2w

class MockAlignedSegment:
    def __init__(self, query_name: str, reference_start: int, query_sequence: str, cigartuples: str):
        self.query_name = query_name
        self.reference_start = reference_start
        self.query_sequence = query_sequence # 0 based
        self.reference_end = reference_start + len(query_sequence)
        self.query_qualities = array("B", [38] * len(query_sequence))
        self.cigartuples = []
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        pysam_cigar = {
            "M": 0,
            "I": 1,
            "D": 2,
            "N": 3,
            "S": 4
            # etc
        }
        for i in list(Cigar(cigartuples).items()):
            self.cigartuples.append((pysam_cigar[i[1]], i[0]))

    def add_indels(self, indels_map):
        cnt = self.reference_start
        for i in self.cigartuples:
            if i[0] == 1: # insert TODO Justify -1
                indels_map.append((self.query_name, self.reference_start, cnt-1, i[1], 0))
            elif i[0] == 2: # del
                indels_map.append((self.query_name, self.reference_start, cnt, 0, 1))
                cnt += i[1]
            else:
                cnt += i[1]

@pytest.mark.parametrize("mArr,spec", [(
    [
        MockAlignedSegment(
            "89.6-2108",
            2291,
            "AAGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATTGAGTTGCCAGGGAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGAGCAGATAGACATAGAAATCTGTGGACATAAAGCTAAAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCA",
            "3M1D61M1I6M1D175M"
        )
    ],[
        "CAGATGATACAGTATTAGAAGAATTGAG-TTGCCAGGGAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGAGCAGATAGACATAGAAATCTGTGGACATAAAGCTAAAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTT"
    ]

)])
def test_some_func(mArr, spec, mocker):

    indels_map = []
    for m in mArr:
        m.add_indels(indels_map)

    mock_samfile = mocker.MagicMock()
    mock_samfile.configure_mock(
        **{
            "fetch.return_value": mArr
        }
    )

    mock_dict = mocker.MagicMock()
    mock_dict.__getitem__.return_value = 42

    window_length = 201
    window_start = 2334 # 0 based

    arr, arr_read_qualities_summary, arr_read_summary, counter = b2w._run_one_window(
        mock_samfile,
        window_start,
        "HXB2-does-not-matter",
        window_length,
        0,
        mock_dict,
        0,
        True,
        indels_map
    )
    print(arr)

    for idx, el in enumerate(arr):
        assert el.split("\n")[1] == spec[idx]