from dataclasses import dataclass

@dataclass(frozen=True)
class ImmutableClass:
    int_number: int
    str_value: str


@dataclass(frozen=False)
class MutableClass:
    int_number: int
    str_value: str


# Testing immutable class
aaa = ImmutableClass(int_number=1, str_value='yo')
bbb = ImmutableClass(int_number=2, str_value='ye')

print(aaa.int_number)
print(bbb.int_number)
print(aaa.str_value)
print(bbb.str_value)

ccc = MutableClass(int_number=1, str_value='yo')
ddd = MutableClass(int_number=2, str_value='ye')

print(ccc.int_number)
print(ddd.int_number)
print(ccc.str_value)
print(ddd.str_value)