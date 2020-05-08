import unittest
import io
import sys
from unittest.mock import Mock, patch
from rnaseq_tools import utils

class EmptyObject:
    pass

class UtilsTester(unittest.TestCase):
    def setUp(self):
        self.patcher = patch('os.path.exists', return_value = True)
        self.patcher = patch('rnaseq_tools.utils.executeSubProcess')
        self.patcher.start()

    def test_getRunNumber(self):
        fastq_path = 'run_1573_s_4_withindex_sequence_GAGTACG.fastq.gz'
        run_number = '1573'
        function_run_number = utils.getRunNumber(fastq_path)
        self.assertEqual(run_number, function_run_number)

    def test_decomposeStatus2Bit(self):
        status = 6
        decomp = [2.0, 1.0]
        function_decomp = utils.decomposeStatus2Bit(status)
        self.assertEqual(decomp, function_decomp)

    def test_makeCombinations(self):
        all_combos = [['a', 'b', 'c'], ['a', 'b'], ['a', 'c'], ['b', 'c']]
        lst = ['a', 'b', 'c']
        function_combos = utils.makeCombinations(lst)
        self.assertEqual(all_combos, function_combos)

    def test_addForwardSlash(self):
        path = '/path/to/file'
        path_with_slash = '/path/to/file/'
        function_path = utils.addForwardSlash(path)
        self.assertEqual(path_with_slash, function_path)

    def test_removeForwardSlash(self):
        path = '/path/to/file/'
        no_end_slash = '/path/to/file'
        func_path = utils.removeForwardSlash(path)
        self.assertEqual(no_end_slash, func_path)

    def test_checkCSV(self):
        file_path = 'path/blah.csv'
        self.assertTrue(utils.checkCSV(file_path))

    def test_checkTSV(self):  # TODO: WRITE THIS. check both return and exception
        pass

    def test_checkExcel(self):  # TODO: WRITE THIS. check both return and exception
        pass

    def test_readInDataframe(self):  # TODO: WRITE THIS. check both return and exception
        pass

    def test_fileBaseName(self):
        file = 'sequence_path.fastq.gz'
        basename = 'sequence_path'
        self.assertEqual(basename, utils.fileBaseName(file))

    def test_pathBaseName(self):
        path = '/path/to/dir/sequence_path.fastq.gz'
        basename = 'sequence_path'
        self.assertEqual(basename, utils.pathBaseName(path))

    def test_dirName(self):
        path_with_file = '/path/to/dir/sequence_path.fastq.gz'
        path_with_filedirname = 'dir'
        self.assertEqual(path_with_filedirname, utils.dirName(path_with_file))
        path_no_file = '/path/to/dir/'
        dir_name = 'dir'
        func_dirname = utils.dirName(path_no_file)
        self.assertEqual(dir_name, func_dirname)

    def test_extractTopmostFiles(self):  # TODO: Write
        pass

    def test_softLinkAndSetAttr(self):
        empty_object = EmptyObject()
        list_dirs = ['dir1', 'dir2']
        origin_dir_path = '/path/from/root'
        soft_link_loc = '/path/to/softlink_here'
        utils.softLinkAndSetAttr(empty_object, list_dirs, origin_dir_path, soft_link_loc)
        self.assertEqual('/path/to/softlink_here/dir1', empty_object.dir1)
        self.assertEqual('/path/to/softlink_here/dir2', empty_object.dir2)

    def test_setAttributes(self):
        new_stdout = io.StringIO()
        sys.stdout = new_stdout

        empty_object = EmptyObject
        empty_object.attr1 = 'first attribute'
        warning_statement = "attr2 not in expected attributes. Not a problem, but also not handled in the objects\n\n"
        utils.setAttributes(empty_object, ['attr1'], {'attr2': 'second attribute'})
        print_statement = new_stdout.getvalue()
        self.assertEqual(empty_object.attr2, 'second attribute')
        self.assertEqual(warning_statement, print_statement)

    def test_userInputCorrectPath(self):
        empty_object = EmptyObject()
        empty_object.attr1 = 'blah'
        utils.userInputCorrectPath('', empty_object, 'attr1', 'new_value')
        self.assertTrue('new_value', empty_object.attr1)

    def test_userInputCorrectAttributeName(self):
        empty_object = EmptyObject()
        empty_object.attr1 = 'blah'
        utils.userInputCorrectAttributeName('correct attr name', empty_object, 'attr1', ['attr_new'])
        self.assertTrue(hasattr(empty_object, 'attr_new') == True)

    def test_createLogger(self):
        logger = utils.createLogger('/dev/null', 'loggertest')
        self.assertLogs(logger, 'WARNING')

    def tearDown(self):
        self.patcher.stop()


if __name__ == '__main__':
    unittest.main()
