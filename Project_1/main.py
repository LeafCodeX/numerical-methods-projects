import os
import pandas as pd
import matplotlib.pyplot as plt
from colorama import Fore, Style


def display_data(data_column, num_entries, project_name):
    for i in range(num_entries):
        if i % 6 == 0 and i != 0:
            print()
        if not pd.isna(data_column.iloc[i]):
            print(f'({i}) ', end="")
            print(Fore.RED + f"{round(data_column.iloc[i], 2)}", end=", ")
            print(Style.RESET_ALL, end="")
        else:
            print(f'({i}) {round(data_column.iloc[i], 2)}', end=", ")
    print()
    separator(project_name)


def separator(project_name):
    alternate_character = "="
    for i in range(len(project_name)):
        print(alternate_character, end="")
        if alternate_character == "=":
            alternate_character = "-"
        else:
            alternate_character = "="
    print()


def exponential_moving_average(data, n, input_column, output_column):
    alpha = 2 / (n + 1)
    for i in range(n, len(data)):
        nominator = 0
        denominator = 0
        for j in range(0, n + 1):
            p = data[input_column].iloc[i - j]
            nominator += p * (1 - alpha) ** j
            denominator += (1 - alpha) ** j
        data.loc[i, f'{output_column}'] = nominator / denominator
    return data


def moving_average_convergence_divergence(data):
    data['MACD'] = data['EMA_12'] - data['EMA_26']
    return data


def intersects_MACD_SIGNAL(data):
    for i in range(1, len(data)):
        if data['MACD'].iloc[i] > data['SIGNAL'].iloc[i] and data['MACD'].iloc[i - 1] < data['SIGNAL'].iloc[i - 1]:
            data.loc[i, 'Buy'] = True
        if data['MACD'].iloc[i] < data['SIGNAL'].iloc[i] and data['MACD'].iloc[i - 1] > data['SIGNAL'].iloc[i - 1]:
            data.loc[i, 'Sell'] = True
    return data


def plot_stock_generator(interval, column1, period, i):
    plt.figure(figsize=(12, 6))
    plt.plot(interval['Date'], interval['Close'], label='Value', color='dodgerblue', linewidth=1)

    plt.title(column1, font='Comic Sans MS', fontweight='bold', fontsize=20)
    plt.xlabel('Date', font='Comic Sans MS', fontsize=14, fontweight='bold')
    plt.ylabel('Price', font='Comic Sans MS', fontsize=14, fontweight='bold')

    plt.xticks(interval['Date'][::period], rotation=45, font='Comic Sans MS', fontweight='bold')
    plt.yticks(font='Comic Sans MS', fontweight='bold')
    plt.tight_layout()
    plt.legend(loc='upper left')
    business = ['Apple_Inc', 'Cd_Project_SA', 'Bitcoin_Cryptocurrency']
    file_name = 'Stock_chart_of_' + str(business[i]) + '_1000_entries.png'
    file_path = os.path.join('Images', file_name)
    plt.savefig(file_path)
    plt.show()


def plot_stock_with_mark_intersections(interval, column1, period, i):
    plt.figure(figsize=(12, 6))
    plt.plot(interval['Date'], interval['Close'], label='Value', color='dodgerblue', linewidth=1)

    plt.scatter(interval['Date'].loc[interval['Buy'] == True], interval['Close'].loc[interval['Buy'] == True],
                color='green', label='Buy', marker='^')
    plt.scatter(interval['Date'].loc[interval['Sell'] == True], interval['Close'].loc[interval['Sell'] == True],
                color='red', label='Sell', marker='v')


    plt.title(column1, font='Comic Sans MS', fontweight='bold', fontsize=20)
    plt.xlabel('Date', font='Comic Sans MS', fontsize=14, fontweight='bold')
    plt.ylabel('Price', font='Comic Sans MS', fontsize=14, fontweight='bold')

    plt.xticks(interval['Date'][::period], rotation=45, font='Comic Sans MS', fontweight='bold')
    plt.yticks(font='Comic Sans MS', fontweight='bold')
    plt.tight_layout()
    plt.legend(loc='upper left')
    business = ['Apple_Inc', 'Cd_Project_SA', 'Bitcoin_Cryptocurrency']
    file_name = str(business[i]) + '_stock_chart_with_mark_intersections_1000_entries.png' if period == 100 else str(
        business[i]) + '_stock_chart_with_mark_intersections_500_entries.png'
    file_path = os.path.join('Images', file_name)
    plt.savefig(file_path)
    plt.show()


def plot_MACD_SIGNAL_buy_sell_generator(interval, title, period, i):
    plt.figure(figsize=(12, 6))
    plt.plot(interval['Date'], interval['MACD'], label='MACD', color='darkgreen', linewidth=1)
    plt.plot(interval['Date'], interval['SIGNAL'], label='SIGNAL', color='orange', linewidth=1)
    plt.scatter(interval['Date'].loc[interval['Buy'] == True], interval['MACD'].loc[interval['Buy'] == True],
                color='green', label='Buy', marker='^')
    plt.scatter(interval['Date'].loc[interval['Sell'] == True], interval['MACD'].loc[interval['Sell'] == True],
                color='red', label='Sell', marker='v')

    plt.title(title, font='Comic Sans MS', fontweight='bold', fontsize=20)
    plt.xlabel('Date', font='Comic Sans MS', fontsize=14, fontweight='bold')
    plt.ylabel('Value', font='Comic Sans MS', fontsize=14, fontweight='bold')

    plt.xticks(interval['Date'][::period], rotation=45, font='Comic Sans MS', fontweight='bold')
    plt.yticks(font='Comic Sans MS', fontweight='bold')
    plt.tight_layout()
    plt.legend(loc='upper left')
    company = ['Apple_Inc', 'Cd_Project_SA', 'Bitcoin_Cryptocurrency']
    file_name = 'MACD_SIGNAL_chart_of_' + str(
        company[i]) + '_1000_entries.png' if period == 100 else 'MACD_SIGNAL_chart_' + str(
        company[i]) + '_500_entries.png'
    file_path = os.path.join('Images', file_name)
    plt.savefig(file_path)
    plt.show()


def algorithm(data, company, advanced):
    start_money = 0.0
    start_actions = 1000
    print(Fore.LIGHTCYAN_EX + 'Company name: ' + str(company) + (
        ' (Simple)' if not advanced else ' (Advanced)') + Style.RESET_ALL)
    print(Fore.LIGHTRED_EX + 'Initial capital: ' + str(start_money) + '$')
    print('Number of initial shares: ' + str(start_actions) + ' => ' + str(start_actions * data['Close'].iloc[0]) + '$')
    print('Number of initial exchange volumes: ' + str(round(start_actions / 1000, 3)) + Fore.RESET)
    if advanced is False:
        start_money, start_actions = simple_algorithm_to_invest_money(data, start_money, start_actions, company, advanced)
    if advanced is True:
        start_money, start_actions = advanced_algorithm_to_invest_money(data, start_money, start_actions, company, advanced)
    print(Fore.LIGHTGREEN_EX + 'Final capital: ' + str(round(start_money, 2)) + '$ => ' + str(
        round(start_money / data['Close'].iloc[-1], 2)) + ' shares')
    print('Number of final shares: ' + str(round(start_actions, 2)) + ' => ' + str(
        round(start_actions * data['Close'].iloc[-1], 2)) + '$')
    print('Number of final exchange volumes: ' + str(round(start_actions / 1000, 3)) + Fore.RESET)


def simple_algorithm_to_invest_money(data, start_money, start_actions, company, advanced):
    money_history_1 = []
    actions_history_1 = []
    dates_1 = []
    for i in range(1, len(data)):
        # MACD cuts from above -> sell actions
        if data['MACD'].iloc[i] >= data['SIGNAL'].iloc[i] and data['MACD'].iloc[i - 1] < data['SIGNAL'].iloc[i - 1]:
            if start_money == 0:
                start_money += start_actions * data['Close'].iloc[i]
                start_actions = 0
            # print("Sell: " + data['Date'].iloc[i])
            # print("Money: " + str(round(start_money, 2)) + "$" + " => " + str(round(start_money / data['Close'].iloc[i], 2)))
            # print("Date: " + data['Date'].iloc[i])
        # MACD cuts from below -> buy action
        elif data['MACD'].iloc[i] <= data['SIGNAL'].iloc[i] and data['MACD'].iloc[i - 1] > data['SIGNAL'].iloc[i - 1]:
            if start_actions == 0:
                start_actions += start_money / data['Close'].iloc[i]
                start_money = 0
            # print("Buy: " + data['Date'].iloc[i])
            # print("Actions: " + str(round(start_actions, 2)) + " => " + str(round(start_actions * data['Close'].iloc[i], 2)) + "$")
            # print("Date: " + data['Date'].iloc[i])
        money_history_1.append(start_money)
        actions_history_1.append(start_actions)
        dates_1.append(data['Date'].iloc[i])
    plot_money_and_actions(dates_1, money_history_1, actions_history_1, company, advanced)
    return start_money, start_actions


def advanced_algorithm_to_invest_money(data, start_money, start_actions, company, advanced):
    money_history_2 = []
    actions_history_2 = []
    dates_2 = []
    for i in range(1, len(data)):
        if data['MACD'].iloc[i - 1] > data['SIGNAL'].iloc[i - 1] and data['MACD'].iloc[i] < data['SIGNAL'].iloc[i] and \
                data['Low'].iloc[i] < data['Low'].iloc[i - 1]:
            if start_money == 0:
                start_money += start_actions * data['Close'].iloc[i]
                start_actions = 0
            # print("Sell: " + data['Date'].iloc[i])
            # print("Money: " + str(round(start_money, 2)) + "$" + " => " + str(round(start_money / data['Close'].iloc[i], 2)))
            # print("Date: " + data['Date'].iloc[i])
        elif data['MACD'].iloc[i - 1] < data['SIGNAL'].iloc[i - 1] and data['MACD'].iloc[i] > data['SIGNAL'].iloc[i] and \
                data['High'].iloc[i] > data['High'].iloc[i - 1]:
            if start_actions == 0:
                start_actions += start_money / data['Close'].iloc[i]
                start_money = 0
            # print("Buy: " + data['Date'].iloc[i])
            # print("Actions: " + str(round(start_actions, 2)) + " => " + str(round(start_actions * data['Close'].iloc[i], 2)) + "$")
            # print("Date: " + data['Date'].iloc[i])
        money_history_2.append(start_money)
        actions_history_2.append(start_actions)
        dates_2.append(data['Date'].iloc[i])
    plot_money_and_actions(dates_2, money_history_2, actions_history_2, company, advanced)
    return start_money, start_actions


def plot_money_and_actions(dates, money_history, actions_history, company, advanced):
    fig, axs = plt.subplots(2, figsize=(12, 12))

    axs[0].plot(dates, actions_history, label='Actions', color='dodgerblue', linewidth=1)
    axs[0].fill_between(dates, actions_history, color='dodgerblue', alpha=0.3)
    axs[0].set_title(str(company) + ' - Shares and money (' + (
        'Simple Algorithm' if not advanced else 'Advanced Algorithm') + ')', font='Comic Sans MS', fontweight='bold',
                     fontsize=20)
    axs[0].set_xlabel('Date', font='Comic Sans MS', fontsize=14, fontweight='bold')
    axs[0].set_ylabel('Amount', font='Comic Sans MS', fontsize=14, fontweight='bold')
    axs[0].legend(loc='upper right')
    axs[0].set_xticks(dates[::240])
    axs[0].set_xticklabels(dates[::240], font='Comic Sans MS', fontweight='bold')
    axs[0].tick_params(axis='both', labelsize=14, labelfontfamily='Comic Sans MS')

    axs[1].plot(dates, money_history, label='Money', color='darkgreen', linewidth=1)
    axs[1].fill_between(dates, money_history, color='darkgreen', alpha=0.3)
    axs[1].set_xlabel('Date', font='Comic Sans MS', fontsize=14, fontweight='bold')
    axs[1].set_ylabel('Value', font='Comic Sans MS', fontsize=14, fontweight='bold')
    axs[1].legend(loc='upper left')
    axs[1].set_xticks(dates[::240])
    axs[1].set_xticklabels(dates[::240], font='Comic Sans MS', fontweight='bold')
    axs[1].tick_params(axis='both', labelsize=14, labelfontfamily='Comic Sans MS')

    business = None
    if company == 'Apple Inc.':
        business = 'Apple_Inc'
    elif company == 'CD Project SA':
        business = 'Cd_Project_SA'
    elif company == 'Bitcoin Cryptocurrency USD':
        business = 'Bitcoin_Cryptocurrency'
    file_name = 'Money_and_actions_of_' + str(business) + ('_simple_algorithm.png' if not advanced else '_advanced_algorithm.png')
    file_path = os.path.join('Images', file_name)
    plt.savefig(file_path)

    plt.tight_layout()
    plt.show()


def main(project_name):
    file_path = ['Data/aapl_us_d.csv', 'Data/cdr_d.csv', 'Data/btc_v_d.csv']
    company_name = ['Apple Inc.', 'CD Project SA', 'Bitcoin Cryptocurrency USD']
    i = 0
    if not os.path.exists('Images'):
        os.makedirs('Images')
    while i < len(file_path):
        file_name = file_path[i]
        data = pd.read_csv(file_name)

        separator(project_name)
        print(Fore.LIGHTBLUE_EX + f'{project_name}' + Style.RESET_ALL)
        separator(project_name)

        # Implementation of the Exponential Moving Average (EMA)
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 1] Calculation of the 12th period EMA:' + Style.RESET_ALL)
        exponential_moving_average(data, 12, 'Close', 'EMA_12')
        display_data(data['EMA_12'], 18, project_name)
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 1] Calculation of the 26th period EMA:' + Style.RESET_ALL)
        exponential_moving_average(data, 26, 'Close', 'EMA_26')
        display_data(data['EMA_26'], 30, project_name)

        # Implementation of the Moving Average Convergence / Divergence indicator
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 2] Calculation of the Moving Average Convergence / Divergence indicator:' +
              Style.RESET_ALL)
        moving_average_convergence_divergence(data)
        display_data(data['MACD'], 30, project_name)

        # Implementation of the SIGNAL
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 3] Calculation of the SIGNAL:' + Style.RESET_ALL)
        exponential_moving_average(data, 9, 'MACD', 'SIGNAL')
        display_data(data['SIGNAL'], 42, project_name)

        # Implementation of the Buy and Sell signals
        intersects_MACD_SIGNAL(data)
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 4] Buy signal:' + Style.RESET_ALL)
        print(data[data['Buy'] == True])
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 4] Sell signal:' + Style.RESET_ALL)
        print(data[data['Sell'] == True])

        # Plotting the stock chart (1000 entries)
        date_start_one = '2020-03-25'
        date_end_one = '2024-03-15'
        interval = data[(data['Date'] >= date_start_one) & (data['Date'] <= date_end_one)]
        plot_stock_generator(interval, str(company_name[i]) + ' stock chart (' + str(date_start_one) + ':' + str(
            date_end_one) + ')', 100, i)
        # Plotting the stock chart with mark intersections (1000 entries)
        plot_stock_with_mark_intersections(interval, str(company_name[i]) + ' stock chart with mark intersections (' + str(
            date_start_one) + ':' + str(date_end_one) + ')', 100, i)
        # Plotting the stock chart with mark intersections (500 entries)
        date_start_two = '2022-03-20'
        date_end_two = '2023-03-15'
        short_interval = data[(data['Date'] >= date_start_two) & (data['Date'] <= date_end_two)]
        plot_stock_with_mark_intersections(short_interval, str(company_name[i]) + ' stock chart with mark intersections (' + str(
            date_start_two) + ':' + str(date_end_two) + ')', 40, i)
        # Plotting the MACD/SIGNAL chart (1000 entries)
        plot_MACD_SIGNAL_buy_sell_generator(interval, 'MACD/SIGNAL chart of ' + str(company_name[i]) + ' (' + str(
            date_start_one) + ':' + str(date_end_one) + ')', 100, i)
        # Plotting the stock chart and the MACD/SIGNAL chart (500 entries)
        plot_MACD_SIGNAL_buy_sell_generator(short_interval, 'MACD/SIGNAL chart of ' + str(company_name[i]) + ' (' + str(
            date_start_two) + ':' + str(date_end_two) + ')', 40, i)

        # Simulation of the stock market actions based on the MACD indicator
        separator(project_name)
        print(Fore.LIGHTBLUE_EX + f'>> [TASK 5] Capital after automatic algorithm based on MACD:' + Style.RESET_ALL)
        advanced = False
        algorithm(interval, company_name[i], advanced)
        advanced = True
        algorithm(interval, company_name[i], advanced)
        separator(project_name)
        i += 1


if __name__ == '__main__':
    main('>> [TITLE] Moving Average Convergence / Divergence (MACD) stock market indicator ')