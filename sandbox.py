from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
import timeit
import cProfile

spe = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/shot133.Lis')
spe.build_spe()
cProfile.run('spe.__build_spe__()', )
t1 = timeit.timeit('spe.__build_spe__()', number=100, globals=globals())
print(t1)

