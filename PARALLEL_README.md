# ECE Forward Modeling - Parallel Computing Guide

이 문서는 ECE forward modeling 코드의 병렬 계산 기능 사용법을 설명합니다.

## 병렬화 개요

기존 코드의 3중 중첩 루프:
```
for cn in range(cnum):           # 채널별 루프
    for i in range(dz.size):     # z 방향 sub-ray
        for j in range(fsub.size): # 주파수별 sub-ray
```

## 병렬화 옵션

### 1. 채널 수준 병렬화 (추천)
- **언제 사용**: 채널 수가 많을 때 (>= 4채널)
- **장점**: 프로세스 간 독립적, 메모리 효율적
- **단점**: 채널 수가 적으면 효과 제한적

```python
# 사용 예시
efm = EceFwdMod()
efm.set_profile(geqdsk_fn, Te_fn, ne_fn)
efm.set_channel(shot, clist)

# 채널별 병렬화 (기본 설정)
results = efm.run(
    parallel=True,
    parallel_mode='channel',  # 채널별 병렬화
    n_processes=4            # 4개 프로세스 사용
)
```

### 2. Sub-ray 수준 병렬화
- **언제 사용**: 채널 수는 적지만 sub-ray가 많을 때 (Nf * Nz가 큰 경우)
- **장점**: 세밀한 병렬화 가능
- **단점**: 메모리 오버헤드 증가

```python
# Sub-ray별 병렬화
results = efm.run(
    parallel=True,
    parallel_mode='sub_ray',  # sub-ray별 병렬화
    n_processes=8            # 8개 스레드 사용
)
```

### 3. 순차 처리 (기존 방식)
```python
# 순차 처리
results = efm.run(parallel=False)
```

## 성능 최적화 팁

### 1. 적절한 프로세스 수 선택
```python
import multiprocessing as mp

# CPU 코어 수 확인
print(f"Available CPU cores: {mp.cpu_count()}")

# 권장 설정
n_processes = min(mp.cpu_count(), len(channel_list))
```

### 2. 병렬화 모드 선택 기준

| 시나리오 | 권장 모드 | 이유 |
|----------|-----------|------|
| 채널 수 > 4, Nf*Nz < 100 | `channel` | 프로세스 생성 비용 대비 효과적 |
| 채널 수 < 4, Nf*Nz > 100 | `sub_ray` | 세밀한 병렬화로 성능 향상 |
| 채널 수 > 8, Nf*Nz > 100 | `channel` | 메모리 효율성 우선 |

### 3. 메모리 사용량 고려
- 채널별 병렬화: 메모리 사용량 = 기본 × 프로세스 수
- Sub-ray별 병렬화: 추가 오버헤드 있음

## 예제 코드

### 기본 사용법
```python
from efm import EceFwdMod

# 모델 초기화
efm = EceFwdMod()
efm.set_profile("data/geqdsk.dat", "data/Te.dat", "data/ne.dat")
efm.set_channel(shot=13728, clist=['ECEI_GT1501', 'ECEI_GT1502'])

# 병렬 실행
Rch, zch, int_meas, tau, rad_temp, abs_temp = efm.run(
    fstart=-0.1, fend=0.1, Nf=5,
    zstart=-10, zend=10, Nz=5,
    parallel=True,
    n_processes=4
)
```

### 성능 비교 테스트
```python
import time

# 순차 처리 시간 측정
start = time.time()
results_seq = efm.run(parallel=False)
seq_time = time.time() - start

# 병렬 처리 시간 측정
start = time.time()
results_par = efm.run(parallel=True, n_processes=4)
par_time = time.time() - start

print(f"Speedup: {seq_time/par_time:.2f}x")
```

## 주의사항

1. **메모리 부족**: 많은 프로세스 사용 시 메모리 부족 가능
2. **I/O 병목**: 파일 읽기가 많으면 병렬화 효과 제한적
3. **디버깅**: 병렬 실행 시 print 출력이 섞일 수 있음
4. **재현성**: 병렬 실행 순서는 비결정적

## 문제 해결

### 메모리 부족 오류
```python
# 프로세스 수를 줄여서 실행
results = efm.run(parallel=True, n_processes=2)
```

### 성능이 오히려 느려지는 경우
```python
# 순차 처리로 변경
results = efm.run(parallel=False)
```

### Import 오류
```python
# parallel_utils.py가 같은 디렉토리에 있는지 확인
# 또는 PYTHONPATH에 추가
```

이 병렬화 구현으로 많은 채널과 sub-ray가 있는 경우 상당한 성능 향상을 기대할 수 있습니다.
