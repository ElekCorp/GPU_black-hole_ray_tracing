project('blackhole', ['cpp', 'cuda'],
default_options : ['cpp_std=c++17','optimization=3','b_ndebug=true'])

add_project_arguments(
  '-m64',
  '-use_fast_math',
  language : 'cuda'
)
add_project_arguments(
  '-m64',
  language: 'cpp'
)


src = ['kernel.cu', 'cli_parser.cpp', 'szinsaver.cpp']
executable('main', src)
