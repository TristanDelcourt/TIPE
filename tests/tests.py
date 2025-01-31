import time
import os
import asyncio
import telegram_send

async def send_message(message):
    await telegram_send.send(messages=[message])


def optimal_b(numbers, bitcount):
    prev = 37_000
    start = time.time()
    os.system(f"../c/bin/factor -b {prev} -s 100000 -q -n {numbers[0]}")
    prev_time = time.time() - start
    
    current = 35_000
    start = time.time()
    os.system(f"../c/bin/factor -b {current} -s 100000 -q -n {numbers[0]}")
    current_time = time.time() - start

    diff = .0005
    while(abs(current_time - prev_time) > diff and abs(current - prev) > 100):
        if(current_time < prev_time):
            if(current < prev):
                current, prev = current - abs(current - prev), current
            else:
                current, prev = (current + prev)//2, current
                
        if(current_time > prev_time):
            if(current > prev):
                current, prev = current + abs(current - prev), current
            else:
                current, prev = (current + prev)//2, current
        
        current = max(current, 50)
        prev_time = current_time
        current_time = 0

        for n in numbers:
            start = time.time()
            os.system(f"../c/bin/factor -b {current} -s 100000 -q -n {n}")
            end = time.time()
            current_time += end - start

        current_time /= len(numbers)
        print(f"{current} | {current_time}")

    return current, current_time

async def main(numbers, bitcount):
    b, t = optimal_b(numbers, bitcount)
    await send_message(f"Tested {len(bit100)} {bitcount} bit numbers, with a optimal B value of {b}, factoring in {t}s")
    print("--------------------------------------------")

    prev = 100_000
    start = time.time()
    os.system(f"../c/bin/factor -b {b} -s {prev} -q -n {numbers[0]}")
    prev_time = time.time() - start
    
    current = 150_000
    start = time.time()
    os.system(f"../c/bin/factor -b {b} -s {current} -q -n {numbers[0]}")
    current_time = time.time() - start

    diff = .0005
    while(abs(current_time - prev_time) > diff and abs(current - prev) > 500):
        if(current_time < prev_time):
            if(current < prev):
                current, prev = current - abs(current - prev), current
            else:
                current, prev = (current + prev)//2, current
                
        if(current_time > prev_time):
            if(current > prev):
                current, prev = current + abs(current - prev), current
            else:
                current, prev = (current + prev)//2, current
        
        current = max(current, 50)
        prev_time = current_time
        current_time = 0

        for n in numbers:
            start = time.time()
            os.system(f"../c/bin/factor -b {b} -s {current} -q -n {n}")
            end = time.time()
            current_time += end - start

        current_time /= len(numbers)
        print(f"{current} | {current_time}")

    await send_message(f"Tested {len(numbers)} {bitcount} bit numbers, with a optimal S value of {current}, factoring in {current_time}s")
    time.sleep(5)
    await send_message(f"*Tested {len(numbers)} {bitcount}, optimal for: b={b}, S={current}*")
    return current, current_time


bit160 = [
    557804934900869097376832247787319559011911034869,
    223239931849808652429153883493331526663330115641,
    173448733838047210325673363532059150117800238227,
    578618697224416436487260319367451944552718797493,
    300028126402977761360766650574333059685903809399
]


bit140 = [
    292580127043266243988988722089282259853953,
    734679417149085432201205638516858819812699,
    113583085715144342210516050744141019877107,
    539390690427396583330171445206318064893331,
    231207090195811462071875779751455825109459
]

# approx. 15200
bit120 = [
    208529505848528486746432688980628747,
    249668092094239699159481381321706631,
    387215876272708326068430335534077811,
    152921157017626342293818827838719171,
    468256643000100204083825032280503379
]

# approx. 5700
bit100 = [
    399163146776588588565277953599,
    744573716629551586149567423413,
    248912655583819057385527989889,
    146852462601173978776099387871,
    289618835801239281914639910817
]


asyncio.run(main(bit140, 140))