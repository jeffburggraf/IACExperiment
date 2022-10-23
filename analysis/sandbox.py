import time
i = 0
lines = []

while True:
    line="\t"*i
    line += f"def f():\n"
    line += "\t"*(i + 1) + "pass"

    lines.append(line)

    out = "\n".join(lines)

    i += 1

    try:
        exec(out)
    except Exception as e:
        print(out)
        print(f"Failed at {len(lines) - 1} nested functions.", e)
        break

    print(len(lines))

    time.sleep(0.3)


